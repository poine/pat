#! /usr/bin/env python

import logging, sys, os, time, math, numpy as np
import rospy, geometry_msgs.msg

import pdb

import pat3.pprz_bag as p3_ppbag
import pat3.algebra as pal
# FIXME should be in ros_pat
import pat3.ros_utils as p3_rpu

# <message name="ROTORCRAFT_FP" id="147">
#       <field name="east"     type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="north"    type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="up"       type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="veast"    type="int32" alt_unit="m/s" alt_unit_coef="0.0000019"/>
#       <field name="vnorth"   type="int32" alt_unit="m/s" alt_unit_coef="0.0000019"/>
#       <field name="vup"      type="int32" alt_unit="m/s" alt_unit_coef="0.0000019"/>
#       <field name="phi"      type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="theta"    type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="psi"      type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="carrot_east"   type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="carrot_north"  type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="carrot_up"     type="int32" alt_unit="m" alt_unit_coef="0.0039063"/>
#       <field name="carrot_psi"    type="int32" alt_unit="deg" alt_unit_coef="0.0139882"/>
#       <field name="thrust"        type="int32"/>
#       <field name="flight_time"   type="uint16" unit="s"/>
# </message>

class MyDataSet(p3_ppbag.DataSet):
    def __init__(self):
        def rec_rotorcraft_fp(msg):
            return [[float(msg.east)*0.0039063, float(msg.north)*0.0039063, float(msg.up)*0.0039063],
                    np.deg2rad([float(msg.phi)*0.0139882, float(msg.theta)*0.0139882, float(msg.psi)*0.0139882])]
        p3_ppbag.DataSet.__init__(self, msg_ids=['ROTORCRAFT_FP'], msg_cbks=[rec_rotorcraft_fp])


class Node(p3_rpu.PeriodicNode):
    def __init__(self, ds):
        self.ds = ds
        self.tf_pub = p3_rpu.TransformPublisher()
        self.marker_pub = p3_rpu.PoseArrayPublisher(dae='quad_murat.dae')
        self.odom_pub = p3_rpu.OdomPublisher()
        self.i=0

    def periodic(self):
        now = rospy.Time.now()
        t_enu, e_ned = self.ds.vals[0][self.i][0], self.ds.vals[0][self.i][1]
        t_ned = [t_enu[1], t_enu[0], -t_enu[2]]
        self.T_b2w = p3_rpu.T_of_t_rpy(t_ned, e_ned)
        #self.T_b2w = pal.T_of_t_eu(t_ned, e_ned)
        self.tf_pub.send_w_enu_to_ned_transform(now)
        self.tf_pub.send_w_ned_to_b_transform(now, self.T_b2w)
        self.marker_pub.publish([self.T_b2w])
        self.odom_pub.publish(self.T_b2w, now)

        self.i += 1
        if self.i >= len(self.ds.vals[0]): self.i=0
        

def main(filename):
    d = MyDataSet()
    d.readPprzBag(filename)
    rospy.init_node('play_rotorcraft_log')
    Node(d).run(25)
    
        
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    filename = '/home/poine/work/pat/data/logs/20_08_04__13_41_08.data'
    main(filename)
