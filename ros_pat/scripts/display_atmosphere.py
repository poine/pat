#! /usr/bin/env python
import numpy as np
import rospy, visualization_msgs.msg, geometry_msgs.msg
import dynamic_reconfigure.server as dyn_rec_srv

import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as pmtf
import pat3.algebra as pal
import pat3.atmosphere as p3_atm

import pat3.ros_utils as p3_rpu

import ros_pat.cfg.atmosphereConfig

import pdb

'''
This is a ROS agent that sends messages for displaying a trajectory in rviz
'''


        
class Agent:

    def __init__(self, atm):
        #self.atm = p3_atm.AtmosphereThermal1()
        #self.atm = p3_atm.AtmosphereThermalMoving()
        self.atm = p3_atm.AtmosphereThermalMulti()
        #self.atm.set_params(xc=15., yc=15., zi=850., wstar=256.)
        #self.atm_pub = p3_rpu.AtmPublisher( self.atm, z0=-14. )
        self.atm_pub = p3_rpu.AtmPointCloudPublisher(self.atm, center=[0, 0, -50], dx=200., dy=200., dz=-100., dens=5.)
        self.tf_pub = p3_rpu.TransformPublisher()
        self.atm_cfg_srv = dyn_rec_srv.Server(ros_pat.cfg.atmosphereConfig,
                                                       self.atm_dyn_cfg_callback)
        
    def atm_dyn_cfg_callback(self, config, level):
        rospy.loginfo(" ATM Reconfigure Request:")
        if level == -1 or level == 1:
            self.atm.set_params(config['thxc'], config['thyc'], config['thzi'], config['thwstar'], config['idx'])
        if  level == -1 or level == 0:
            xc = [config['xc'], config['yc'], -config['hc']]
            print xc, config['dx'], config['dy'], -config['dh'], config['spacing']
            self.atm_pub.update_grid(xc, config['dx'], config['dy'], -config['dh'], config['spacing'])
        return config
     
    def periodic(self):
        self.atm_pub.update_atm(self.atm, rospy.Time.now().to_sec())
        self.atm_pub.publish(self.atm)
        self.tf_pub.publish(rospy.Time.now(), T_w2b=None)
    
    def run(self):
        rate = rospy.Rate(2.)
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass

def main():
    rospy.init_node('display_atmosphere')
    #traj_name = rospy.get_param('~traj_name', 'circle4')
    #if rospy.get_param('~list_available', False):
    #    rospy.loginfo('  available trajectories\n{}'.format(pmtf.list()))
    #try:
        #rospy.loginfo('  Loading trajectory: {}'.format(traj_name))
    atm = None
        #rospy.loginfo('  Description: {}'.format(desc))
    Agent(atm).run()
    #except KeyError:
    #    rospy.loginfo('  Unkwown trajectory: {}'.format(traj_name))
    #    rospy.loginfo( '  available trajectories\n{}'.format(pmtf.list()))
        
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
