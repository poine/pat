#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

''' Plot of aerodynmics '''
import pat3.utils as p3_u

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.plot_utils as p3_pu
import pat3.frames as p3_fr


# redo with equilibrium Usfc
def plot_polar(dm):
    alpha1, alpha2 = np.deg2rad([-1., 15])
    beta, rvel, Usfc = 0, [0, 0,0], [0, 0, 0, 0]
    alphas = np.linspace(alpha1, alpha2, 100)
    ACoefs = np.array([p1_fw_dyn.get_f_aero_coef(alpha, beta, rvel, Usfc, dm.P) for alpha in alphas])
    CLs, CDs = ACoefs[:,0], ACoefs[:,2] 
    #plt.plot(alphas, CLs)
    LovDs = CLs/CDs; i_max = np.argmax(LovDs)
    print('max L over D  {:.1f} (alpha {:.2f} deg)'.format(LovDs[i_max], np.rad2deg(alphas[i_max])))
    plt.plot(CDs, CLs)
    p3_pu.decorate(plt.gca(), xlab='CD', ylab='CL')
    plt.plot



    
    

def trim_throttle(dm):
    trim_args={'h':0, 'va':12, 'gamma':0}
    Xe, Ue = dm.trim(trim_args, debug=True)
    trim_args={'h':0, 'va':12, 'throttle':0.1762}
    Xe, Ue = dm.trim(trim_args, debug=True)
    trim_args={'h':0, 'va':12, 'throttle':0.}
    Xe, Ue = dm.trim(trim_args, debug=True)
#    pdb.set_trace()
    
def plot_trims(dm, throttle=0., force_recompute=False, nvs=10, nhs=10):
    vs = np.linspace(5, 30, nvs, endpoint=True)
    hs = np.linspace(0, 10000, nhs, endpoint=True)
    trims = np.zeros((len(hs), len(vs), 2))
    alphas, gammas, zdots = np.zeros((len(hs), len(vs))), np.zeros((len(hs), len(vs))), np.zeros((len(hs), len(vs))) 
    filename = '/tmp/foo2.npz'
    if force_recompute or not os.path.exists(filename):
        for i,v in enumerate(vs):
            for j,h in enumerate(hs):
                Xe, Ue = dm.trim({'va': vs[i], 'h': hs[i], 'throttle': throttle})
                alphas[j, i], gammas[j, i] = Xe[dm.sv_alpha], Xe[dm.sv_theta]-Xe[dm.sv_alpha]
                Xe_ee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(Xe)
                zdots[j, i] = Xe_ee[ p3_fr.SixDOFEuclidianEuler.sv_zd]
        np.savez(filename, alphas=alphas, gammas=gammas, zdots=zdots)
    else:
        data =  np.load(filename)
        alphas, gammas, zdots =  data['alphas'], data['gammas'], data['zdots']

    ax = plt.subplot(1,3,1)
    alphas_deg = np.rad2deg(alphas)
    v = np.linspace(np.min(alphas_deg), np.min([np.max(alphas_deg), 15]), 32, endpoint=True)
    plt.contour(vs, hs, alphas_deg, v, linewidths=0.5, colors='k')
    plt.contourf(vs, hs, alphas_deg, v, cmap=plt.cm.jet)
    plt.colorbar(ticks=v)
    p3_pu.decorate(ax, title='alpha', xlab='velocity (m/s)', ylab='height (m)')

    ax = plt.subplot(1,3,2)
    gammas_deg = np.rad2deg(gammas)
    v = np.linspace(np.min(gammas_deg), np.max(gammas_deg), 32, endpoint=True)
    plt.contour(vs, hs, gammas_deg, v, linewidths=0.5, colors='k')
    plt.contourf(vs, hs, gammas_deg, v, cmap=plt.cm.jet)
    plt.colorbar(ticks=v)
    p3_pu.decorate(ax, title='gamma', xlab='velocity (m/s)', ylab='height (m)')

    ax = plt.subplot(1,3,3)
    v = np.linspace(np.min(zdots), 1., 32)#np.max(zdots), 32, endpoint=True)
    plt.contour(vs, hs, zdots, v, linewidths=0.5, colors='k')
    plt.contourf(vs, hs, zdots, v, cmap=plt.cm.jet)
    plt.colorbar(ticks=v)
    p3_pu.decorate(ax, title='zdot', xlab='velocity (m/s)', ylab='height (m)')
    
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/skywalker_x8.xml')
    dm = p1_fw_dyn.DynamicModel(param_filename)
    plot_polar(dm)
    #trim_throttle(dm)
    plot_trims(dm)
    plt.show()

        
if __name__ == "__main__":
    main()

