import rclpy
import rclpy.node, rclpy.clock
import tf_transformations
import tf2_ros
from tf2_ros.static_transform_broadcaster import StaticTransformBroadcaster
from geometry_msgs.msg import TransformStamped
import visualization_msgs.msg, geometry_msgs.msg

import numpy as np

import pdb

class TransformPublisher:
    def __init__(self, parent_node):
        self.tf_broadcaster = tf2_ros.TransformBroadcaster(parent_node)
        self.tf_static_broadcaster = StaticTransformBroadcaster(parent_node)
        #self.send_w_enu_to_ned_transform(rclpy.clock.Clock().now())

    def send_w_enu_to_ned_transform(self, t):
        R_enu2ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        T_enu2ned = np.eye(4); T_enu2ned[:3,:3] = R_enu2ned
        self.send_transform("w_enu", "w_ned", t, T_enu2ned, static=True)

    def send_w_ned_to_b_transform(self, t, T_w2b):
        self.send_transform("w_ned", "b_frd", t, T_w2b)
        
    def send_transform(self, f1, f2, t, T_f1tof2, static=False):
        tf_msg = TransformStamped()
        tf_msg.header.frame_id = f1
        tf_msg.child_frame_id = f2
        tf_msg.header.stamp = t.to_msg()
        _r = tf_msg.transform.rotation
        _r.x, _r.y, _r.z, _r.w = tf_transformations.quaternion_from_matrix(T_f1tof2)
        _t = tf_msg.transform.translation
        _t.x, _t.y, _t.z = T_f1tof2[:3,3]
        if static:
            self.tf_static_broadcaster.sendTransform(tf_msg)
        else:
            self.tf_broadcaster.sendTransform(tf_msg)


        
class OdomPublisher:
    def __init__(self, parent_node, topic='ros_pat/odom', frame='w_ned'):
        pass

def _position_and_orientation_from_T(p, q, T):
    p.x, p.y, p.z = T[:3, 3]
    q.x, q.y, q.z, q.w = tf_transformations.quaternion_from_matrix(T)    

def T_of_t_q(t, q):
    T = tf_transformations.quaternion_matrix(q)
    T[:3,3] = t
    return T


class MarkerArrayPublisher:
    def __init__(self, parent_node, topic, mtype=2,
                 meshes=['package://ros2_pat/meshes/quadcopter.dae'],
                 colors=[[0.2, 1., 0.2, 0.5]], scales=[(1., 1., 1.)], frame_id="w_ned"):
        self.marker_pub = parent_node.create_publisher(visualization_msgs.msg.MarkerArray, topic, 2)
        self.msg = visualization_msgs.msg.MarkerArray()
        for i, (mesh, color, scale) in enumerate(zip(meshes, colors, scales)):
            marker = visualization_msgs.msg.Marker()
            marker.header.frame_id = frame_id
            marker.type = mtype #marker.MESH_RESOURCE
            marker.action = marker.ADD
            marker.id = i
            #marker.text = "{}".format(i)
            marker.scale.x, marker.scale.y, marker.scale.z = scale
            marker.color.r, marker.color.g, marker.color.b, marker.color.a  = color
            #p = marker.pose.position; p.x, p.y, p.z = 0.3, 0, 0
            marker.mesh_resource = mesh
            marker.mesh_use_embedded_materials = True
            self.msg.markers.append(marker)

    def publish(self, T_ned2bs, delete=False):
        for marker, T_ned2b in zip(self.msg.markers, T_ned2bs):
            marker.action = marker.DELETE if delete else marker.ADD
            _position_and_orientation_from_T(marker.pose.position, marker.pose.orientation, T_ned2b)
        #pdb.set_trace()
        self.marker_pub.publish(self.msg)


class MarkerPublisher:
    def __init__(self, parent_node, topic, mtype=2,
                 mesh='package://ros2_pat/meshes/quadcopter.dae',
                 color=[0.2, 1., 0.2, 0.5], scale=(1., 1., 1.), frame_id="w_ned"):
         self.marker_pub = parent_node.create_publisher(visualization_msgs.msg.Marker, topic, 2)
         self.msg = visualization_msgs.msg.Marker()
         self.msg.header.frame_id = frame_id
         self.msg.type = mtype #marker.MESH_RESOURCE
         self.msg.action = self.msg.ADD
         self.msg.scale.x, self.msg.scale.y, self.msg.scale.z = scale
         self.msg.color.r, self.msg.color.g, self.msg.color.b, self.msg.color.a  = color
         self.msg.mesh_resource = mesh
         #self.msg.mesh_use_embedded_materials = True
         
    def publish(self, stamp, T_ned2b, delete=False):
        self.msg.header.stamp = stamp.to_msg()
        _position_and_orientation_from_T(self.msg.pose.position, self.msg.pose.orientation, T_ned2b)
        self.marker_pub.publish(self.msg)
        

class TrajPublisher:
    def __init__(self, traj, parent_node, topic, frame_id="w_ned", color=[0.2, 1., 0.2, 0.5]):
        self.marker_pub = parent_node.create_publisher(visualization_msgs.msg.Marker, topic, 2)
        self.msg = visualization_msgs.msg.Marker()
        self.msg.header.frame_id = frame_id
        self.msg.type = 4 # line strip FIXME, get the constant
        self.msg.action = self.msg.ADD
        self.msg.id = 1 
        self.msg.pose.orientation.w = 1.
        self.msg.scale.x, self.msg.scale.y, self.msg.scale.z = 0.01, 0.01, 0.01
        self.msg.color.r, self.msg.color.g, self.msg.color.b, self.msg.color.a = color
        for l in np.linspace(0,1, 100):
        #for _p in [[0.,0.,0.], [0.,1.,0.], [2.,1.,0.]]:
            #pdb.set_trace()
            p = geometry_msgs.msg.Point()
            p.x, p.y, p.z = traj._g.get(l)[:3,0]
            self.msg.points.append(p)
        
    def publish(self):
        #print(self.msg, '\n')
        self.marker_pub.publish(self.msg)
