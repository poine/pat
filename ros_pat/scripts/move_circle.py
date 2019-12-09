#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg, rosbag
import dynamic_reconfigure.client
import matplotlib.pyplot as plt

import pat3.ros_utils as p3_rpu
import pat3.plot_utils as p3_pu


def plot_bag(bag_path='2019-12-03-13-53-31.bag'):
    bag = rosbag.Bag(bag_path, "r")
    _t, _z, _zd = [], [], []
    for topic, msg, ts in bag.read_messages(topics=['/guidance/status']):
        _t.append(ts.to_sec())
        _z.append(msg.h)
        _zd.append(msg.hdot)

    _t = np.array(_t)
    _z = np.array(_z)
    _zd = np.array(_zd)
    ax = plt.subplot(2,1,1)
    plt.plot(_t, _z)
    p3_pu.decorate(ax, title='height', xlab='time in s', ylab='z m')
    ax = plt.subplot(2,1,2)
    plt.plot(_t, _zd)
    p3_pu.decorate(ax, title='vz', xlab='time in s', ylab='vz in m/s', legend=True)
    plt.show()




class Agent(p3_rpu.PeriodicNode):
    def __init__(self):
        p3_rpu.PeriodicNode.__init__(self, 'move_circle')
        self.c0 = np.array([0, 0, 0])       # ENU ?
        self.v =  np.array([0., 0.1, 0])   # 
        self.guid_dyn_rec_client = dynamic_reconfigure.client.Client("/real_time_sim", timeout=1, config_callback=self.guid_config_callback)
        self.t0 = rospy.Time.now().to_sec()
        
    def guid_config_callback(self, config):
        #print config
        pass
        
    def periodic(self):
        #print 'hello'
        dt = rospy.Time.now().to_sec() - self.t0
        circle_c = self.c0 + dt*self.v
        self.guid_dyn_rec_client.update_configuration({"circle_xc":circle_c[0], "circle_yc":circle_c[1]})
        

def main():
    #Agent().run(freq=5.)
    plot_bag()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
