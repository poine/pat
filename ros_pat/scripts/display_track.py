#! /usr/bin/env python
import os, numpy as np
import rospy, rospkg


import pat3.ros_utils as pru

class Agent:
    def __init__(self, track_path):
        self.track_pub = pru.TrackPublisher(track_path)

    def periodic(self):
        #now = rospy.Time.now()
        self.track_pub.publish()

    def run(self):
        rate = rospy.Rate(2.)
        t0 = rospy.Time.now().to_sec()
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass  
        

        

def main():
    rospy.init_node('display_track')
    track_name = rospy.get_param('~track_name', 'track1')
    #ppwd = rospkg.RosPack().get_path('ros_pat')
    #track_path = os.path.join(ppwd, 'data/tracks/{}.yaml'.format(track_name))
    track_path = '/home/poine/work/pat/data/tracks/track1.yaml'
    rospy.loginfo('  Loading track: {}'.format(track_name))
    Agent(track_path).run()
    

if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
