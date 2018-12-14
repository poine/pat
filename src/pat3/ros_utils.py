import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg, visualization_msgs.msg, yaml
import pat3.algebra as pal
import pdb

def _position_and_orientation_from_T(p, q, T):
    p.x, p.y, p.z = T[:3, 3]
    q.x, q.y, q.z, q.w = tf.transformations.quaternion_from_matrix(T)


class TransformPublisher:
    def __init__(self):
        self.tfBcaster = tf2_ros.TransformBroadcaster()
        self.tfBuffer  = tf2_ros.Buffer()
        self.tfLstener = tf2_ros.TransformListener(self.tfBuffer)

    def publish(self, t, T_w2b):
        self.send_w_enu_to_ned_transform(t)
        self.send_w_ned_to_b_transform(t, T_w2b) 
    
    def send_w_enu_to_ned_transform(self, t):
        R_enu2ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        T_enu2ned = np.eye(4); T_enu2ned[:3,:3] = R_enu2ned
        self.send_transform("w_enu", "w_ned", t, T_enu2ned)

    def send_w_ned_to_b_transform(self, t, T_w2b):
        self.send_transform("w_ned", "b_frd", t, T_w2b)

    def send_transform(self, f1, f2, t, T_f1tof2):
        tf_msg = geometry_msgs.msg.TransformStamped()
        tf_msg.header.frame_id = f1
        tf_msg.child_frame_id = f2
        tf_msg.header.stamp = t
        _r = tf_msg.transform.rotation
        _r.x, _r.y, _r.z, _r.w = tf.transformations.quaternion_from_matrix(T_f1tof2)
        _t = tf_msg.transform.translation
        _t.x, _t.y, _t.z = T_f1tof2[:3,3]
        self.tfBcaster.sendTransform(tf_msg)


        
class MarkerArrayPublisher:
    def __init__(self, topic, meshes, colors=[[0.2, 1., 0.2, 0.5]]):
        self.meshes = meshes
        self.pub = rospy.Publisher(topic, visualization_msgs.msg.MarkerArray, queue_size=1)
        self.msg = visualization_msgs.msg.MarkerArray()
        for i, (mesh, color) in enumerate(zip(meshes, colors)):
            marker = visualization_msgs.msg.Marker()
            marker.header.frame_id = "w_ned"
            marker.type = marker.MESH_RESOURCE
            marker.action = marker.ADD
            marker.id = i
            marker.text = "{}".format(i)
            marker.scale.x, marker.scale.y, marker.scale.z = 1., 1., 1.
            marker.color.r, marker.color.g, marker.color.b, marker.color.a  = color
            marker.mesh_resource = mesh
            marker.mesh_use_embedded_materials = True
            self.msg.markers.append(marker)
        
    def publish(self, T_ned2bs):
        for marker, T_ned2b in zip(self.msg.markers, T_ned2bs):
            _position_and_orientation_from_T(marker.pose.position, marker.pose.orientation, T_ned2b)
        self.pub.publish(self.msg)

class PoseArrayPublisher(MarkerArrayPublisher):
    def __init__(self):
        MarkerArrayPublisher.__init__(self, '/pat/vehicle_marker',  ["package://ros_pat/media/quad.dae"])


class QuadAndRefPublisher(MarkerArrayPublisher):
    def __init__(self):
       meshes = ["package://ros_pat/media/quad.dae", "package://ros_pat/media/quad.dae"]
       colors = [[0.2, 1., 0.2, 0.75], [0.7, 0.7, 0.7, 0.5]]
       MarkerArrayPublisher.__init__(self, '/pat/vehicle_marker',  meshes, colors)


class TrackPublisher(MarkerArrayPublisher):
    def __init__(self, path):
        with open(path, 'r') as stream:
            d = yaml.load(stream)
        self.poses, colors, meshes = [], [], []
        for name, obj in d.iteritems():
            pos = np.array(obj['position'])
            ori = np.array(obj['orientation'])
            self.poses.append(pal.T_of_t_eu(pos, ori))
            meshes.append(obj['mesh'])
            colors.append(obj['color'])
        MarkerArrayPublisher.__init__(self, '/pat/track_poles', meshes, colors)

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
        #marker.points.append(marker.points[0])

        
    def publish(self):
        self.traj_pub.publish(self.traj_msg)
    
