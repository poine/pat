#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
Learning a linear pitch controller
 see: http://recherche.enac.fr/~drouin/autom_avancee/tuto_ann_2/
'''

import os, os.path, sys, pickle, logging
import numpy as np, matplotlib.pyplot as plt
import keras
import control
import pdb

LOG = logging.getLogger(__name__)

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.utils as p3_u
import pat3.frames as p3_fr

import so_lti
import utils as ut

def make_or_load_plant_id_training_set(_ac, dm, Xe, Ue, dt=0.005, nb_samples=10000, force_recompute=False):
    filename = '/tmp/pat_plant_training_set_{}.pkl'.format(_ac)
    if force_recompute or not os.path.exists(filename):
        theta_e, q_e, ele_e = Xe[p3_fr.SixDOFAeroEuler.sv_theta], Xe[p3_fr.SixDOFAeroEuler.sv_q], Ue[dm.iv_de()]
        LOG.info('creating plant-id training dataset')
        
        delta_vas    = np.random.uniform(-0.5, 0.5, nb_samples)#np.zeros(nb_samples)
        delta_alphas = np.random.uniform(np.deg2rad(-0.5),  np.deg2rad(0.5), nb_samples)#np.zeros(nb_samples)
        delta_thetas = np.random.uniform(np.deg2rad(-20),  np.deg2rad(20),  nb_samples)
        delta_qs     = np.random.uniform(np.deg2rad(-200), np.deg2rad(200), nb_samples)
        delta_deles  = np.random.uniform(np.deg2rad(-20),  np.deg2rad(20),  nb_samples)
        Xkp1s = []
        for d_vak, d_alphak, d_thetak, d_qk, d_delek in zip(delta_vas, delta_alphas, delta_thetas, delta_qs, delta_deles):
            Xk = np.array(Xe) + [0, 0, 0,   d_vak, d_alphak, 0,   0, d_thetak, 0,   0, d_qk, 0]
            #Uk = np.array(Ue) + [0, 0, delta_delek, 0]#, 0]
            Uk = np.array(Ue); Uk[dm.iv_de()]+= d_delek
            dm.reset(Xk, 0, Uk)
            Xkp1 = dm.run(dt=dt, tf=dt, U=Uk, atm=None)
            Xkp1s.append(Xkp1)
        Xkp1s = np.array(Xkp1s)
        _input = np.vstack([delta_thetas, delta_qs, delta_deles]).T
        _output = np.vstack([Xkp1s[:,p3_fr.SixDOFAeroEuler.sv_theta]-theta_e, Xkp1s[:,p3_fr.SixDOFAeroEuler.sv_q]-q_e]).T
        with open(filename, "wb") as f:
            pickle.dump([_input, _output], f)
            LOG.info('saved to {}'.format(filename))

    else:
        LOG.info('Loading plant-id training dataset from {}'.format(filename))
        with open(filename, "rb") as f:
            _input, _output = pickle.load(f)
    return _input, _output

    
def identify_plant(_ac, v_trim=12., dt=0.005, epochs=100, force_retrain=False, force_remake_training_set=False ):
    param_filename = p3_u.pat_ressource('data/vehicles/{}.xml'.format(_ac))
    plant_ann_filename='/tmp/pat_pil_plant_ann_{}.h5'.format(_ac)
    dm = p1_fw_dyn.DynamicModel(param_filename)
    trim_args = {'h':0, 'va':v_trim, 'gamma':0.}
    Xe, Ue = dm.trim(trim_args, report=True)
    if force_retrain or force_remake_training_set or not os.path.exists(plant_ann_filename):
        _input, _output = make_or_load_plant_id_training_set(_ac, dm, Xe, Ue, force_recompute=force_remake_training_set)

        # Build the plant identification ANN
        plant_i = keras.layers.Input((3,), name ="plant_i") # x1_k, x2_k, u_k
        plant_l = keras.layers.Dense(2, activation='linear', kernel_initializer='uniform',
                                     input_shape=(3,), use_bias=False, name="plant")
        plant_o = plant_l(plant_i)
        plant_ann = keras.models.Model(inputs=plant_i, outputs=plant_o)
        plant_ann.compile(loss='mean_squared_error', optimizer='adam')
        plant_ann.fit(_input, _output, epochs=epochs, batch_size=32,  verbose=1, shuffle=True)
        plant_ann.save(plant_ann_filename)
    else:
         plant_ann = keras.models.load_model(plant_ann_filename)
    return dm, Xe, Ue, plant_ann

def report_plant_id(plant, Xe, Ue, plant_ann, dt=0.005):
    Ac, Bc = plant.get_jacobian(Xe, Ue)
    s = p3_fr.SixDOFAeroEuler
    A1c = np.array([[Ac[s.sv_theta, s.sv_theta], Ac[s.sv_theta, s.sv_q]],
                    [Ac[s.sv_q, s.sv_theta], Ac[s.sv_q, s.sv_q]]])
    B1c = np.array([[Bc[s.sv_theta, plant.iv_de()]],
                    [Bc[s.sv_q,     plant.iv_de()]]])
    ct_ss = control.ss(A1c, B1c, [[1, 0]], [[0]])
    dt_ss = control.sample_system(ct_ss, dt, method='zoh') #  ‘matched’, ‘tustin’, ‘zoh’
    fmt_arg = {'precision':6, 'suppress_small':True}
    LOG.info('Real plant jacobian:\n{}\n{}'.format(np.array2string(dt_ss.A, **fmt_arg), np.array2string(dt_ss.B)))
    LOG.info(' modes: {}'.format(np.linalg.eig(dt_ss.A)[0]))
    w_identified = plant_ann.get_layer(name="plant").get_weights()[0]
    #print('Identified weights: \n{}'.format(w_identified))
    A_id = np.matrix(w_identified[:2,:].T, dtype=np.float64)
    B_id = np.matrix(w_identified[2:,:].T, dtype=np.float64)
    LOG.info('Identified plant:\n{}\n{}'.format(np.array2string(A_id, **fmt_arg), np.array2string(B_id, **fmt_arg)))
    LOG.info(' modes: {}'.format(np.linalg.eig(A_id)[0]))
  

def build_ctl_plus_plant_ann(plant_ann):
    # build ann for controller
    ref_input = keras.layers.Input((4,), name ="ctl_ref_i")   # Xr1_kp1, Xr2_kp1, Xr1_k, Xr2_k         
    state_input = keras.layers.Input((2,), name ="ctl_x_i")   # X1_k, X2_k
    ctl_input = keras.layers.concatenate([ref_input, state_input])
    ctl_l = keras.layers.Dense(1, activation='linear', kernel_initializer='uniform',
                                   input_shape=(6,), use_bias=False, name="ctl")
    ctl_output = ctl_l(ctl_input)
    # build ann for plant
    plant_input =  keras.layers.concatenate([state_input, ctl_output]) # X1_k, X2_k, U_k
    plant_l = plant_ann.get_layer(name="plant")
    plant_l.trainable = False
    plant_output = plant_l(plant_input)
    # build the (controller + plant) network
    full_ann = keras.models.Model(inputs=[ref_input, state_input], outputs=plant_output)
    full_ann.compile(loss='mean_squared_error', optimizer='adam')#, metrics=['accuracy'])
    #full_ann.summary()
    return full_ann 
    
def train_controller(dt, plant_ann, omega_ref=5., xi_ref=0.9, omega_err=10., xi_err=0.7, epochs=20,
                     force_retrain=False, full_ann_filename='/tmp/pat_pil_full_ann.h5'):
    ''' Train the controller '''
    if force_retrain or not os.path.isfile(full_ann_filename): 
        full_ann = build_ctl_plus_plant_ann(plant_ann)
        # build LTI systems for reference model and tracking error dynamics
        ref = so_lti.CCPlant(omega_ref, xi_ref, dt=dt)
        track_err = so_lti.CCPlant(omega_err, xi_err, dt=dt)
        # simulate and record reference model trajectory
        time, Xr, Ur, desc = so_lti.make_or_load_training_set(ref, ut.CtlNone(), True)
        _len = len(Xr)-1
        # simulate and record tracking errors
        eps_k = np.random.uniform(low=-1., high=1., size=(_len,2)) # tracking error
        eps_kp1 = np.array([np.dot(track_err.Ad, _eps_k) for _eps_k in eps_k])
        # format training data
        _input = [np.zeros((_len, 4)), np.zeros((_len, 2))] 
        _output = np.zeros((_len, 2))
        for k in range(_len):
            _input[0][k] = [Xr[k+1, 0], Xr[k+1, 1], Xr[k, 0], Xr[k, 1]]
            _input[1][k] = Xr[k] + eps_k[k]
            _output[k] = Xr[k+1] + eps_kp1[k]
        full_ann.fit(_input, _output, epochs=epochs, batch_size=16, verbose=1, shuffle=True)
        full_ann.save(full_ann_filename)
    else:
        # load a previously trained ANN
        full_ann = keras.models.load_model(full_ann_filename)
        ref = so_lti.CCPlant(omega_ref, xi_ref, dt=dt)
    # Build the controller-only network
    ctl_i = keras.layers.Input((6,), name ="ctl_i") # Xr1_kp1, Xr2_kp1, Xr1_k, Xr2_k , X1_k, X2_k
    ctl_o = full_ann.get_layer(name="ctl")(ctl_i)
    ann_ctl = keras.models.Model(inputs=ctl_i, outputs=ctl_o)
    #ann_ctl.summary()
    #w_identified = full_ann.get_layer(name="ctl").get_weights()
    #print('Identified ctl: \n{}'.format(w_identified))
    return ref, ann_ctl


def report_controller(ann_ctl, filename='/tmp/ctl.yaml'):
    k_identified = ann_ctl.get_layer(name="ctl").get_weights()[0]
    htkp1, hqkp1, htk, hqk, k_theta, k_q = k_identified.squeeze() 
    H, K = [htkp1, hqkp1, htk, hqk], [-k_theta, -k_q] # FIXME, why - ???
    K.append(-0.005) # integrator....
    txt = (  "H_pitch:\n"
           + "  rows: 1\n"
           + "  cols: 4\n"
           + "  data: [" + ", ".join(["{:.12f}".format(i) for i in H]) + "]\n"
           + "K_pitch:\n"
           + "  rows: 1\n"
           + "  cols: 3\n"
           + "  data: [" + ", ".join(["{:.12f}".format(i) for i in K]) + "]\n")
    with open(filename, 'w') as f:
        f.write(txt)
    LOG.info('saved controller to {}'.format(filename))
    print('{}'.format(txt))

def controller_compute(ref, ann_ctl, Yc_k, X_k, Xr_k):
    ''' compute reference model and controller outputs '''
    Xr_kp1 = ref.disc_dyn(Xr_k, Yc_k)
    Uk = ann_ctl.predict(np.array([[Xr_kp1[0], Xr_kp1[1], Xr_k[0], Xr_k[1], X_k[0], X_k[1]]]))
    return Xr_kp1, Uk

def test_controller(ref, ann_ctl, plant, Xe, Ue, plant_ann):
    ''' run the controller on the real world plant and plot results '''
    time =  np.arange(0., 30.01, plant.dt)
    sp, Xr = ut.step_vec(time, a0=np.deg2rad(-1), a1=np.deg2rad(1), dt=20., t0=5.), np.zeros((len(time),2))
    #sp, Xr = np.array([p3_u.step(_t, a=ampl, p=_p, dt=_p/4.) for _t in time]), np.zeros((len(time),2))
    X, U = np.zeros((len(time),2)), np.zeros((len(time),1))
    _s = p3_fr.SixDOFAeroEuler
    Xp, Up = np.zeros((len(time), _s.sv_size)), np.zeros((len(time), plant.input_nb()))
    theta_e, q_e = Xe[_s.sv_theta], Xe[_s.sv_q]
    X[0] = [0.2, 0]
    Xp[0] = np.array(Xe); Xp[0, _s.sv_theta] += np.deg2rad(2)
    Xr[0,0] = Xp[0, _s.sv_theta]-theta_e
    plant.reset(Xp[0])
    #lin_plant = so_lti.Plant(np.eye(2), np.zeros((2,1)), 0.01)
    #_w = plant_ann.get_layer(name="plant").get_weights()[0]
    #lin_plant.Ad = _w[:2,:].T
    #lin_plant.Bd = _w[2:,:].T
    ki, sum_theta_err = 0.0075, 0
    for k in range(0,len(time)-1):
        #Xr[k+1], U[k] = controller_compute(ref, ann_ctl, [sp[k]], X[k], Xr[k])
        Xr[k+1], U[k] = controller_compute(ref, ann_ctl, [sp[k]], [Xp[k,_s.sv_theta]-theta_e, Xp[k, _s.sv_q]], Xr[k])
        sum_theta_err += Xp[k,_s.sv_theta]-theta_e-sp[k]
        Up[k] = Ue; Up[k,plant.iv_de()] += U[k]+ki*sum_theta_err
        #X[k+1] = lin_plant.disc_dyn(X[k], U[k])
        Xp[k+1] = plant.run(plant.dt, time[k+1], Up[k], atm=None)
    Up[-1] = Up[-2]
    #so_lti.plot(time, X, U, Xr)
    plant.plot_trajectory(time, Xp, Up, label='airplane')
    plt.subplot(5,3,8)
    plt.plot(time, np.rad2deg(theta_e+sp), label='setpoint')
    plt.plot(time, np.rad2deg(theta_e+Xr[:,0]), label='reference')
    plt.legend()
    plt.subplot(5,3,11)
    plt.plot(time, np.rad2deg(q_e+Xr[:,1]), label='reference')
    plt.legend()
        
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    force_plant_training_set = '-fpts' in sys.argv
    force_plant_id = '-fpid' in sys.argv
    force_ctl_id = '-fcid' in sys.argv
    #_ac, v_trim, dt = 'stormr_avl', 15., 0.005
    _ac, v_trim, dt = 'cularis', 15., 0.005
    # Identify plant
    plant, Xe, Ue, plant_ann = identify_plant(_ac, v_trim, dt, 100, force_plant_id, force_plant_training_set)
    report_plant_id(plant, Xe, Ue, plant_ann, dt)
    # Train the controller on the identified plant
    ref, ann_ctl = train_controller(dt, plant_ann, force_retrain=force_ctl_id, epochs=100)
    report_controller(ann_ctl)
    test_controller(ref, ann_ctl, plant, Xe, Ue, plant_ann)
    plt.show()
    
    #pdb.set_trace()
