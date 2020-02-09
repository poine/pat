#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

''' Fitting netto vario parameters '''
import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.atmosphere as p3_atm
import pat3.plot_utils as p3_pu
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn


def test_phi0(atm, dm, h=0):
    vas = np.arange(6, 20, 1.)
    vzs = []
    for va in vas:
        Xae_e, Uae_e = dm.trim({'h':h, 'va':va, 'throttle':0}, report=False, debug=False)
        vzs.append(-Xae_e[p3_fr.SixDOFEuclidianEuler.sv_zd])

    rho = p3_atm.get_rho(h)
    K = 2.*dm.P.m*dm.P.g/rho/dm.P.Sref
    H = np.array([[va**3/K, K/va] for va in vas])
    CD0, B = np.dot(np.linalg.pinv(H), vzs)
    print('K {} CD0 {} B {}'.format(K, CD0, B))

    vas2 = np.arange(6, 20, 0.1)
    vzs2 = [va**3/K*CD0+K/va*B for va in vas2]
    
    plt.plot(vas, vzs, '*')
    plt.plot(vas2, vzs2)
    p3_pu.decorate(plt.gca(), xlab='va (m/s)', ylab='vz (m/s)')
    plt.show()

def test_phi(atm, dm):
    phis = np.deg2rad(np.arange(-30, 31, 1.))
    h, va = 0, 14

    Xes, Ues = [], []
    for phi in phis:
        Xe, Ue = dm.foo_trim_aero_banked_cst_throttle(h=h, va=va, throttle=0., flaps=0., phi=phi, debug=True)
        Xes.append(Xe); Ues.append(Ue)

    Xes, Ues = np.array(Xes), np.array(Ues)
    #pdb.set_trace()
    # trims = np.array([dm.trim_aero_banked_cst_throttle(h=h, va=va, throttle=0., flaps=0., phi=phi, debug=True)
    #          for phi in phis])
    # trims_Xe = []
    # for (gamma, elevator, aileron, rudder, alpha), phi in zip(trims, phis):
    #     theta = gamma+alpha
    #     psid = dm.P.g*np.tan(phi)/va
    #     p, q, r = 0., np.sin(phi)*np.cos(theta)*psid, np.cos(phi)*np.cos(theta)*psid
    #     trims_Xe.append(dm.from_six_dof_aero_euler([0., 0., -h, va, alpha, 0., phi, theta, 0., p, q, r]))
    # trims_Xe = np.array(trims_Xe)
        
    plt.figure()
    plt.subplot(4,1,1)
    #plt.plot(np.rad2deg(phis), np.rad2deg(trims[:,0]), label='gamma')
    #p3_pu.decorate(plt.gca(), title='gamma', xlab='phi (deg)', ylab='deg')
    plt.plot(np.rad2deg(phis), np.rad2deg(Xes[:,p3_fr.SixDOFEuclidianEuler.sv_theta]), label='theta')
    p3_pu.decorate(plt.gca(), title='theta', xlab='phi (deg)', ylab='deg')
    plt.subplot(4,1,2)
    plt.plot(np.rad2deg(phis), np.rad2deg(Ues[:,dm.iv_de()]), label='ele')
    p3_pu.decorate(plt.gca(), title='ele', xlab='phi (deg)', ylab='deg')
    plt.subplot(4,1,3)
    plt.plot(np.rad2deg(phis), np.rad2deg(Ues[:,dm.iv_da()]), label='ail')
    plt.plot(np.rad2deg(phis), np.rad2deg(Ues[:,dm.iv_dr()]), label='rud')
    p3_pu.decorate(plt.gca(), title='ail/rud', xlab='phi (deg)', ylab='deg', legend=True)
    plt.subplot(4,1,4)
    #plt.plot(np.rad2deg(phis), np.rad2deg(trims[:,4]), label='alpha')
    #p3_pu.decorate(plt.gca(), title='alpha', xlab='phi (deg)', ylab='deg')
    plt.plot(np.rad2deg(phis), Xes[:,p3_fr.SixDOFEuclidianEuler.sv_zd], label='zd')
    p3_pu.decorate(plt.gca(), title='zd', xlab='phi (deg)', ylab='m/s')
    plt.show()

def test_phi2(atm, dm, compute=False): 
    savefile_name = '/tmp/pat_glider_trims.npz'
    h = 0.
    if compute:
        phis = np.deg2rad(np.arange(-45, 46, 1.)) 
        vas  = np.arange(6, 20, 1.)
        vzs = np.zeros((len(vas), len(phis)))
        for iphi, phi in enumerate(phis):
            for iva, va in enumerate(vas):
                Xe, Ue = dm.foo_trim_aero_banked_cst_throttle(h=h, va=va, throttle=0., flaps=0., phi=phi, debug=True)
                vzs[iva, iphi] =  Xe[p3_fr.SixDOFEuclidianEuler.sv_zd]

        np.savez(savefile_name, phis=phis, vas=vas, vzs=vzs)
        print('saved {}'.format(savefile_name))
    else:
        _data =  np.load(savefile_name)
        phis, vas, vzs = [_data[k] for k in ['phis', 'vas', 'vzs']]
        
    # fit
    rho = p3_atm.get_rho(h)
    K = 2.*dm.P.m*dm.P.g/rho/dm.P.Sref
    H,Y = [],[]
    for iphi, phi in enumerate(phis):
        for iva, va in enumerate(vas):
            CL = K/va**2
            #H = np.array([[va**3/K, K/va] for va in vas])
            H.append([va/CL, va*CL/np.cos(phi)**2])
            Y.append(vzs[iva, iphi])
    H=np.array(H)
    Y=np.array(Y)
    CD0, B = np.dot(np.linalg.pinv(H), Y)
    #pdb.set_trace()
    print('K {} CD0 {} B {}'.format(K, CD0, B))
    
    plt.subplot(1,2,1)
    for iva, va in enumerate(vas):
        plt.plot(np.rad2deg(phis), vzs[iva,:], label='va {}'.format(va))
    p3_pu.decorate(plt.gca(), title='zd', xlab='phi (deg)', ylab='m/s', legend=True)
    plt.subplot(1,2,2)
    for iphi, phi in enumerate(phis):
        plt.plot(vas, vzs[:,iphi], label='phi {}'.format(np.rad2deg(phi)))
    p3_pu.decorate(plt.gca(), title='zd', xlab='va (m/s)', ylab='m/s', legend=True)

    vzs2 = np.zeros((len(vas), len(phis)))
    for iphi, phi in enumerate(phis):
        for iva, va in enumerate(vas):
            CL = K/va**2
            vzs2[iva, iphi] = va*(CD0/CL+B*CL/np.cos(phi)**2)
            
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    _phis, _vas = np.meshgrid(phis, vas)
    surf = ax.plot_surface(_vas, np.rad2deg(phis), vzs)
    surf2 = ax.plot_surface(_vas, np.rad2deg(phis), vzs2)
    ax.set_xlabel('va (m/s)')
    ax.set_ylabel('phi (deg)')
    ax.set_zlabel('vz (m/s)')
    
    
    plt.show()
     
def main(h=0):
    atm = p3_atm.AtmosphereCstWind([0., 0., 0.])
    dm = p1_fw_dyn.DynamicModel_ee()
    #test_phi0(atm, dm)
    #test_phi(atm, dm)
    test_phi2(atm, dm, compute=False)
    
    
    

        
if __name__ == "__main__":
    main()
