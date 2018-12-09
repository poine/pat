import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg, visualization_msgs.msg

import pdb

def _position_and_orientation_from_T(p, q, T):
    p.x, p.y, p.z = T[:3, 3]
    q.x, q.y, q.z, q.w = tf.transformations.quaternion_from_matrix(T)


    
class MarkerArrayPublisher:
    def __init__(self, topic, meshes):
        self.meshes = meshes
        self.pub = rospy.Publisher(topic, visualization_msgs.msg.MarkerArray, queue_size=1)
        self.msg = visualization_msgs.msg.MarkerArray()
        for i, mesh in enumerate(meshes):
            marker = visualization_msgs.msg.Marker()
            marker.header.frame_id = "w_ned"
            marker.type = marker.MESH_RESOURCE
            marker.action = marker.ADD
            marker.id = i
            marker.text = "{}".format(i)
            marker.scale.x, marker.scale.y, marker.scale.z = 1., 1., 1.
            marker.color.r, marker.color.g, marker.color.b, marker.color.a  = 0.2, 1., 0.2, 0.5
            marker.mesh_resource = mesh
            marker.mesh_use_embedded_materials = True
            self.msg.markers.append(marker)
            
        
    def publish(self, T_ned2bs):
        for marker, T_ned2b in zip(self.msg.markers, T_ned2bs):
            _position_and_orientation_from_T(marker.pose.position, marker.pose.orientation, T_ned2b)
        self.pub.publish(self.msg)

class PoseArrayPublisher(MarkerArrayPublisher):
    def __init__(self):
        MarkerArrayPublisher.__init__(self, '/pat/reference_marker',  ["package://smocap/meshes/quad2.dae"])

class TrackPublisher(MarkerArrayPublisher):
    def __init__(self):
        self.poses = [np.eye(4), np.eye(4)]
        self.poses[0][0, 3] = -1.
        self.poses[1][0, 3] =  1.
        meshes = ["package://smocap/meshes/fpv_pole.dae", "package://smocap/meshes/fpv_pole.dae"]
        MarkerArrayPublisher.__init__(self, '/pat/track_poles', meshes)

    def publish(self):
        MarkerArrayPublisher.publish(self, self.poses)

class TrajectoryPublisher:

    def __init__(self, traj):
        self.traj_pub = rospy.Publisher('/pat/trajectory', visualization_msgs.msg.MarkerArray, queue_size=1)
        self.traj_msg = visualization_msgs.msg.MarkerArray()
        t0, t1, dt = 0., traj.duration, 0.05
        time = np.arange(t0, t1, dt)
        marker = visualization_msgs.msg.Marker()
        marker.header.frame_id = "w_ned"
        marker.type = marker.LINE_STRIP
        marker.action = marker.ADD
        marker.id = 0
        marker.text = 'trajectory'
        marker.scale.x, marker.scale.y, marker.scale.z = 0.01, 0.01, 0.01
        marker.color.r, marker.color.g, marker.color.b, marker.color.a  = 1, 0, 0, 0.5
        marker.pose.orientation.w = 1.0
        marker.pose.position.x = 0
        marker.pose.position.y = 0
        marker.pose.position.z = 0
        self.traj_msg.markers.append(marker)
        for t in time:
            p = geometry_msgs.msg.Point()
            p.x, p.y, p.z = traj.get(t)[:3, 0]
            marker.points.append(p)
        marker.points.append(marker.points[0])

        
    def publish(self):
        self.traj_pub.publish(self.traj_msg)
    
