#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg
import dynamic_reconfigure.server as dyn_rec_srv
import dynamic_reconfigure.client

import pdb

import pat3.algebra as pal
import pat3.utils as pmu
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d
import pat3.atmosphere as p3_atm
import pat3.ros_utils as p3_rpu
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.guidance as p3_guid

import ros_pat.cfg.guidanceConfig
import ros_pat.cfg.atmosphereConfig
import ros_pat.msg

class Agent:
    def __init__(self, _traj=None, time_factor=1.):
        traj = p3_traj3d.CircleRefTraj(c=[0, 0, -20], r=15)
        #traj = p3_fw_guid.LineRefTraj(p1=[0, 0, 0], p2=[50, 50, 0])
        #traj = p3_fw_guid.SquareRefTraj()
        #traj = p3_fw_guid.OctogonRefTraj()
        self.time_factor = time_factor
        param_filename='/home/poine/work/pat/data/vehicles/cularis.xml'
        _fdm = p1_fw_dyn.DynamicModel(param_filename)
        _fms = p3_guid.FMS(_fdm, {'h':0, 'va':10, 'gamma':0})
        #_ctl = p3_fw_guid.Guidance(_fdm, traj, {'h':0, 'va':10, 'gamma':0})
        #_ctl = p3_guid.GuidanceThermal(_fdm, traj, {'h':0, 'va':10, 'gamma':0})
        #_atm = p3_atm.AtmosphereCstWind([0, 0, -1])
        #_atm = p3_atm.AtmosphereThermal1()
        _atm = p3_atm.AtmosphereThermalMulti()
        #_atm = p3_atm.AtmosphereRidge()

        
        self.sim = pmu.Sim(_fdm, _fms, _atm)
        self.tf_pub = p3_rpu.TransformPublisher()
        self.carrot_pub = p3_rpu.MarkerPublisher('guidance/goal', 'w_ned', scale=(0.25, 0.25, 0.5))
        self.traj_pub = p3_rpu.TrajectoryPublisher2(traj, ms=0.2)
        #self.pose_pub = p3_rpu.PoseArrayPublisher(dae='glider.dae')
        self.pose_pub = p3_rpu.PoseArrayPublisher(dae='ds_glider_full.dae')
        self.atm_pub = None#p3_rpu.AtmPointCloudPublisher(_atm)
        self.status_pub = p3_rpu.GuidanceStatusPublisher()
        self.atm_disp_dyn_rec_client = dynamic_reconfigure.client.Client("display_atmosphere", timeout=1, config_callback=self.atm_config_callback)
        #self.my_own_reconfigure_client = dynamic_reconfigure.client.Client("real_time_sim", timeout=1, config_callback=None)

        self.cur_thermal_id = 0
        self.dyn_cfg_srv = dyn_rec_srv.Server(ros_pat.cfg.guidanceConfig,
                                                   self.dyn_cfg_callback)
        # we can not have two servers... fuck...
        #self.atm_cfg_srv = dyn_rec_srv.Server(ros_pat.cfg.atmosphereConfig,
        #                                      self.atm_dyn_cfg_callback)#, subname='foo')
        self.nav_goal_sub = rospy.Subscriber("/move_base_simple/goal", geometry_msgs.msg.PoseStamped, self.nav_goal_callback)
        self.joy_sub = rospy.Subscriber('/joy', sensor_msgs.msg.Joy, self.joy_callback)
        self.periodic_cnt = 0

    def joy_callback(self, msg):
        #self.last_input = rospy.get_rostime()
        phi_sp, theta_sp = -msg.axes[2], msg.axes[1]
        self.sim.ctl.guidances[self.sim.ctl.mod_auto1].set_setpoints(phi_sp, theta_sp)
        
    def nav_goal_callback(self, msg):
        # frame id = enu
        x, y = round(msg.pose.position.x), round(msg.pose.position.y)
        z = round(float(self.sim.fdm.X[p1_fw_dyn.sv_z]))
        print('nav_goal_callback {} {} {}'.format(x, y, z))
        self.dyn_cfg_srv.update_configuration({'circle_xc': x, 'circle_yc': y, 'circle_zc': z})
        self.dyn_cfg_srv.update_configuration({'fms_mode':p3_guid.FMS.mod_circle})
        
    def atm_config_callback(self, config):
        #rospy.loginfo("Config set to {int_param}, {double_param}, {str_param}, {bool_param}, {size}".format(**config))
        #print config
        print('atm_config_callback')
        
    # def atm_dyn_cfg_callback(self, config, level):
    #     rospy.loginfo(" ATM Reconfigure Request:")
    #     self.sim.atm.set_params(config['thxc'], config['thyc'], config['thzi'], config['thwstar'])
    #     return config

    
    def dyn_cfg_callback(self, config, level):
        rospy.loginfo(" Guidance Reconfigure Request: level {}".format(level))
        self.sim.ctl.set_va_sp(config['vsp'])
      
        self.sim.ctl.set_mode(config['fms_mode'])
        _g = self.sim.ctl.guidance()
        if self.sim.ctl.mode == self.sim.ctl.mod_centering or self.sim.ctl.mode == self.sim.ctl.mod_circle:
            self.traj_pub.update_ref_traj(_g.traj)
        #self.traj_pub.update_ref_traj(self.sim.ctl.traj)
        if level == 1 or level == -1:
            self.dyn_cfg_atmosphere(config)
        if level == 2 or level == -1:
            self.dyn_cfg_guid_circle(config)
        if level == 3 or level == -1:
            self.dyn_cfg_guid_centering(config)
        return config 

    def dyn_cfg_guid_circle(self, config, set_current_z=True):
        print('dyn_cfg_guid_circle'.format(**config))
        _e, _n, _u = config['circle_xc'], config['circle_yc'], config['circle_zc']
        if set_current_z:
            _u = config['circle_zc'] = float(self.sim.fdm.X[p1_fw_dyn.sv_z])
            #_u = self.sim.fdm.X[p1_fw_dyn.sv_z]
        self.sim.ctl.guidances[self.sim.ctl.mod_circle].set_params((_n, _e, _u), config['circle_r'])
        if self.traj_pub is not None:
            self.traj_pub.update_ref_traj(self.sim.ctl.guidance().traj)

    def dyn_cfg_guid_centering(self, config):
        self.sim.ctl.guidances[self.sim.ctl.mod_centering].set_param(config['climbing_gain'])
        self.sim.ctl.guidances[self.sim.ctl.mod_centering].set_radius(config['climbing_radius'])
            
    def dyn_cfg_atmosphere(self, config):
        print('dyn_cfg_atmosphere {idx}'.format(**config))
        if config['idx'] == self.cur_thermal_id:
            print 'params'
            self.sim.atm.set_params(config['thxc'], config['thyc'], config['thzi'], config['thwstar'], config['idx'])
        else:
            print 'idx'
            self.cur_thermal_id = config['idx']
            params = self.sim.atm.get_params(config['idx'])
            # config does not want np.float64
            config['thxc'], config['thyc'], config['thzi'], config['thwstar'] = float(params[0]), float(params[1]), params[2], params[3] 
        if self.atm_pub is not None: self.atm_pub.update_atm(self.sim.atm)
        self.atm_disp_dyn_rec_client.update_configuration({"thxc":config['thxc'], "thyc":config['thyc'], "thzi":config['thzi'], "thwstar":config['thwstar'], "idx":config['idx']})


    
    def periodic(self):
        now = rospy.Time.now()
        t_sim = (now.to_sec() - self.t0)*self.time_factor
        self.sim.run(t_sim)
        self.tf_pub.publish(now, self.sim.fdm.T_w2b)
        #self.pose_pub.publish([self.sim.fdm.T_w2b, self.sim.ctl.T_w2b_ref])
        self.pose_pub.publish([self.sim.fdm.T_w2b])
        self.carrot_pub.publish(self.sim.ctl.carrot())
        if self.traj_pub is not None and len(self.sim.Xs) > 0:# and self.periodic_cnt % 4 == 1:
            #cval = np.array(self.sim.Xs)[:, p3_fr.SixDOFAeroEuler.sv_z]
            cval = -np.array(self.sim.Xees)[:, p3_fr.SixDOFEuclidianEuler.sv_zd]
            if self.sim.ctl.mode == self.sim.ctl.mod_centering:
                self.traj_pub.update_ref_traj(self.sim.ctl.guidance().traj)
            self.traj_pub.publish(self.sim.Xs, self.sim.Xees, cval)
        if self.atm_pub is not None and self.periodic_cnt % 12 == 0: self.atm_pub.publish(self.sim.atm)
        self.status_pub.publish(self)
        self.periodic_cnt += 1
        
    def run(self):
        rate = rospy.Rate(20.)
        self.t0 = rospy.Time.now().to_sec()
        self.sim.reset(0, self.sim.ctl.Xe)
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass  

def main():
    rospy.init_node('real_time_sim')
    #traj_name = rospy.get_param('~traj_name', 'c42')
    #if rospy.get_param('~list_available', False):
    #    rospy.loginfo('  available trajectories\n{}'.format(pmtf.list()))
    #rospy.loginfo('  Loading trajectory: {}'.format(traj_name))
    #traj, desc = pmtf.get(traj_name)
    time_factor = rospy.get_param('~time_factor', 1.)
    rospy.loginfo('  Running simulation with time factor: {}'.format(time_factor))
    Agent(None, time_factor).run()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
