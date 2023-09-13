#
# This is some outdated support for ROS
# This is kept for reference and will be removed shortly
# Support for ROS2 is now in ros2_pat directory
#

import numpy as np
import matplotlib.pyplot as plt
import rospy, tf2_ros, tf
import std_msgs.msg, geometry_msgs.msg, sensor_msgs.point_cloud2, sensor_msgs.msg, visualization_msgs.msg, sensor_msgs.msg, nav_msgs.msg
import yaml
import pat3.algebra as pal
import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.atmosphere as p3_atm
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import ros_pat.msg
import pdb

#
# ROS utilities. Send and receive our messages
#


def _position_and_orientation_from_T(p, q, T):
    p.x, p.y, p.z = T[:3, 3]
    q.x, q.y, q.z, q.w = tf.transformations.quaternion_from_matrix(T)

# Transforms
def T_of_t_rpy(t, rpy):
    T = tf.transformations.euler_matrix(rpy[0], rpy[1], rpy[2], 'sxyz')
    T[:3,3] = t
    return T

def t_rpy_of_T(T):
    t = T[:3,3]
    rpy = tf.transformations.euler_from_matrix(T, 'sxyz')
    return t, rpy

def T_of_t_q(t, q):
    T = tf.transformations.quaternion_matrix(q)
    T[:3,3] = t
    return T

def T_of_t_R(t, R):
    T = np.eye(4); T[:3,:3] = R; T[:3,3] = t
    return T

def T_of_t_r(t, r):
    R, _ = cv2.Rodrigues(r)
    return T_of_t_R(t, R)

def tq_of_T(T):
    return T[:3, 3], tf.transformations.quaternion_from_matrix(T)

def tR_of_T(T):
    return T[:3,3], T[:3,:3]

def tr_of_T(T):
    ''' return translation and rodrigues angles from a 4x4 transform matrix '''
    r, _ = cv2.Rodrigues(T[:3,:3])
    return T[:3,3], r.squeeze()


def transform(a_to_b_T, p_a):
    return np.dot(a_to_b_T[:3,:3], p_a) + a_to_b_T[:3,3]


class TransformPublisher:
    def __init__(self):
        self.tfBcaster = tf2_ros.TransformBroadcaster()
        # FIXME, use that for ENU to NED
        self.static_tfBcaster = tf2_ros.StaticTransformBroadcaster()
        self.tfBuffer  = tf2_ros.Buffer()
        self.tfLstener = tf2_ros.TransformListener(self.tfBuffer)

    def publish(self, t, T_w2b):
        self.send_w_enu_to_ned_transform(t)
        if T_w2b is not None: self.send_w_ned_to_b_transform(t, T_w2b) 
    
    def send_w_enu_to_ned_transform(self, t):
        R_enu2ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        T_enu2ned = np.eye(4); T_enu2ned[:3,:3] = R_enu2ned
        self.send_transform("w_enu", "w_ned", t, T_enu2ned, static=True)

    def send_w_ned_to_b_transform(self, t, T_w2b):
        self.send_transform("w_ned", "b_frd", t, T_w2b)

    def send_w_enu_to_b_transform(self, t, T_w2b, bf='b_frd'):
        self.send_transform("w_enu", bf, t, T_w2b)

    def send_b_to_a_transform(self, t, T_b2a):
        self.send_transform("b_frd", "a_ab", t, T_b2a)
        
    def send_transform(self, f1, f2, t, T_f1tof2, static=False):
        tf_msg = geometry_msgs.msg.TransformStamped()
        tf_msg.header.frame_id = f1
        tf_msg.child_frame_id = f2
        tf_msg.header.stamp = t
        _r = tf_msg.transform.rotation
        _r.x, _r.y, _r.z, _r.w = tf.transformations.quaternion_from_matrix(T_f1tof2)
        _t = tf_msg.transform.translation
        _t.x, _t.y, _t.z = T_f1tof2[:3,3]
        if static:
            self.static_tfBcaster.sendTransform(tf_msg)
        else:
            self.tfBcaster.sendTransform(tf_msg)


        
class MarkerArrayPublisher:
    def __init__(self, topic, meshes, colors=[[0.2, 1., 0.2, 0.5]], scales=[(1., 1., 1.)], frame_id="w_ned"):
        self.meshes = meshes
        self.pub = rospy.Publisher(topic, visualization_msgs.msg.MarkerArray, queue_size=1)
        self.msg = visualization_msgs.msg.MarkerArray()
        for i, (mesh, color, scale) in enumerate(zip(meshes, colors, scales)):
            marker = visualization_msgs.msg.Marker()
            marker.header.frame_id = frame_id
            marker.type = marker.MESH_RESOURCE
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
        self.pub.publish(self.msg)

class PoseArrayPublisher(MarkerArrayPublisher):
    def __init__(self, topic='/pat/vehicle_marker', dae='quad.dae', scales=[(1., 1., 1.)], frame_id="w_ned"):
        MarkerArrayPublisher.__init__(self, topic,  ["package://ros_pat/media/{}".format(dae)], scales=scales, frame_id=frame_id)


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


# 4D traj...
class TrajectoryPublisher:

    def __init__(self, traj, ms=0.01):
        self.traj_pub = rospy.Publisher('/pat/trajectory', visualization_msgs.msg.MarkerArray, queue_size=1)
        self.traj_msg = visualization_msgs.msg.MarkerArray()
        t0, t1, dt = 0., traj.duration, 0.05
        time = np.arange(t0, t1, dt)
        marker = visualization_msgs.msg.Marker()
        marker.header.frame_id = "w_ned"
        marker.type = marker.LINE_STRIP
        marker.action = marker.ADD
        marker.id = 0
        marker.scale.x = ms
        marker.color.r, marker.color.g, marker.color.b, marker.color.a  = 1, 0, 0, 0.5
        marker.pose.orientation.w = 1.0
        p = marker.pose.position; p.x, p.y, p.z = 0, 0, 0
        self.traj_msg.markers.append(marker)
        for t in time:
            p = geometry_msgs.msg.Point()
            p.x, p.y, p.z = traj.get(t)[:3, 0]
            marker.points.append(p)
        #marker.points.append(marker.points[0])# close?

        
    def publish(self):
        self.traj_pub.publish(self.traj_msg)

#
# 3D traj
#
class TrajectoryPublisher2:

    def __init__(self, ref_traj, ms=0.01):
        self.traj_pub = rospy.Publisher('/pat/trajectory', visualization_msgs.msg.MarkerArray, queue_size=1)
        self.traj_msg = visualization_msgs.msg.MarkerArray()
        # reference trajectory
        self.marker_ref = self._create_marker(0, ms, rgba=(1., 1., 0, 0.5))
        self.traj_msg.markers.append(self.marker_ref)
        # reference ground track
        self.marker_ref_gt = self._create_marker(1, ms, rgba=(0.8, 0.8, 0, 0.5))
        self.traj_msg.markers.append(self.marker_ref_gt)
        self.update_ref_traj(ref_traj)
        # vehicle trajectory
        self.marker_v = self._create_marker(2, ms/2, rgba=(1, 1, 0, 0.5))
        self.traj_msg.markers.append(self.marker_v)

    def update_ref_traj(self, ref_traj):
        self.marker_ref.points, self.marker_ref_gt.points = [], []
        for _p in ref_traj.get_points():
            self.marker_ref.points.append(geometry_msgs.msg.Point(*_p))
            p_ref_gc = geometry_msgs.msg.Point(*_p); p_ref_gc.z = 0
            self.marker_ref_gt.points.append(p_ref_gc)
        
    def _create_marker(self, _id, _ms, rgba=(1, 0, 0, 0.5)):
        m = visualization_msgs.msg.Marker()
        m.header.frame_id = "w_ned"
        m.type, m.action, m.id = m.LINE_STRIP, m.ADD, _id
        m.scale.x = _ms
        m.color.r, m.color.g, m.color.b, m.color.a  = rgba
        o = m.pose.orientation; o.w, o.x, o.y, o.z = 1.0, 0, 0, 0
        return m

    def publish(self, Xv=None, Xvee=None, cval=None):
        # publish vehicle trajectory
        if Xv is not None:
            self.marker_v.points = []
            if len(Xv) > 0:
                for X in Xv:
                    p = geometry_msgs.msg.Point(*X[p3_fr.SixDOFAeroEuler.sv_slice_pos])
                    self.marker_v.points.append(p)
                if cval is not None:
                    cmap_name =  ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                                  'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'][10]
                    cval = (cval - cval.min())/cval.ptp()   # Normalize
                    colors = plt.get_cmap(cmap_name)(cval)
                    self.marker_v.colors = [std_msgs.msg.ColorRGBA(*c) for c in  colors]
                
            
        self.traj_pub.publish(self.traj_msg)
    


class MarkerPublisher:
    def __init__(self, topic, ref_frame, scale=(0.05, 0.05, 0.5), argb=(1., 0., 1., 0.)):
        self.carrot_pub = rospy.Publisher(topic, visualization_msgs.msg.Marker, queue_size=1)
        
        self.carrot_msg = visualization_msgs.msg.Marker()
        self.carrot_msg.header.frame_id=ref_frame
        self.carrot_msg.type = visualization_msgs.msg.Marker.CYLINDER
        p = self.carrot_msg.pose.position; p.x, p.y, p.z = 0, 0, 0.025
        o = self.carrot_msg.pose.orientation; o.x, o.y, o.z, o.w = 0, 0, 0, 1
        s = self.carrot_msg.scale; s.x, s.y, s.z = scale
        c = self.carrot_msg.color; c.a, c.r, c.g, c.b = argb

    def publish(self, carrot_pos):
        self.carrot_msg.header.stamp = rospy.Time.now()
        p = self.carrot_msg.pose.position; p.x, p.y, p.z = carrot_pos[0], carrot_pos[1],carrot_pos[2]
        self.carrot_pub.publish(self.carrot_msg)

#
# Atmosphere
#

class AtmPublisher:
    def __init__(self, atm, topic="pat/atmosphere", z0=None):
       
        self.pub = rospy.Publisher(topic, visualization_msgs.msg.MarkerArray, queue_size=1)
        self.msg = visualization_msgs.msg.MarkerArray()
        if z0 is None:
            x, y, z, wx, wy, wz = self.sample(atm)
        else:
            x, y, z, wx, wy, wz = self.sample(atm, x_r=np.arange(-60., 60., 2.5), y_r=np.arange(-60., 60., 2.5),
                                              z_r=np.array([z0]))
        
        w = np.vstack((wx[np.newaxis],wy[np.newaxis],wz[np.newaxis]))
        w_norm = np.linalg.norm(w, axis = 0)
        w_normalized = (w_norm.ravel() - w_norm.min()) / w_norm.ptp()
        cmap_name =  ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                      'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'][10]
        colors = plt.get_cmap(cmap_name)(w_normalized)
        nx, ny, nz = x.shape
        idx=0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    w = np.array([wx[ix, iy, iz], wy[ix, iy, iz], wz[ix, iy, iz]])
                    if np.linalg.norm(w) > 0.15:
                        marker = visualization_msgs.msg.Marker()
                        marker.id = idx
                        marker.header.frame_id = "w_ned"
                        #marker.type = marker.ARROW
                        marker.type = marker.SPHERE
                        marker.action = marker.ADD
                        s = marker.scale; s.x, s.y, s.z = 0.5, 0.5, 0.5
                        c = marker.color; c.r, c.g, c.b, c.a = colors[idx]
                        marker.pose.orientation.w = 1.0
                        p = marker.pose.position; p.x, p.y, p.z = x[ix, iy, iz], y[ix, iy, iz], z[ix, iy, iz]
                        #marker.points=[geometry_msgs.msg.Point(0, 0, 0),
                        #               geometry_msgs.msg.Point(*w)]
                    
                        self.msg.markers.append(marker)
                    idx += 1

    def sample(self, atm, x_r=np.arange(-60., 60., 5), y_r=np.arange(-60., 60., 5), z_r=np.arange(0., -150, -10)):
        x, y, z = np.meshgrid(x_r, y_r, z_r)
        wx, wy, wz = np.meshgrid(x_r, y_r, z_r)
        nx, ny, nz = x.shape
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    pos = [x[ix, iy, iz], y[ix, iy, iz], z[ix, iy, iz]]
                    wx[ix, iy, iz], wy[ix, iy, iz], wz[ix, iy, iz] = atm.get_wind(pos, t=0)
        return x, y, z, wx, wy, wz
                    
                    
    def publish(self):
        self.pub.publish(self.msg)

        

class SimplePublisher(rospy.Publisher):
    def __init__(self, topic, msg_class, what, qs=1, latch=False):
        rospy.loginfo(' {} publishing on {}'.format(what, topic))
        rospy.Publisher.__init__(self, topic, msg_class, queue_size=qs, latch=latch)
        self.msg_class = msg_class

# https://github.com/spillai/pybot/blob/master/pybot/externals/ros/pointclouds.py        
class AtmPointCloudPublisher(SimplePublisher):
    def __init__(self, atm, topic="pat/atmosphere_pc", center=[0, 0, -75], dx=100., dy=100., dz=-150., dens=5. ):
        SimplePublisher.__init__(self, topic, sensor_msgs.msg.PointCloud2, 'simulator') 
        self.header = std_msgs.msg.Header()
        self.header.frame_id = "w_ned"
        self.msg_fields = [
            sensor_msgs.msg.PointField('x', 0, sensor_msgs.msg.PointField.FLOAT32, 1),
            sensor_msgs.msg.PointField('y', 4, sensor_msgs.msg.PointField.FLOAT32, 1),
            sensor_msgs.msg.PointField('z', 8, sensor_msgs.msg.PointField.FLOAT32, 1),
            sensor_msgs.msg.PointField('intensity', 12, sensor_msgs.msg.PointField.FLOAT32, 1)]
        self.update_grid(center, dx, dy, dz, dens)
        self.update_atm(atm)

    def update_grid(self, center, dx, dy, dz, dens):
        print('update grid')
        xc, yc, zc =  center
        self.x_range = np.arange(xc-dx/2, xc+dx/2,  dens)
        self.y_range = np.arange(yc-dy/2, yc+dy/2,  dens)
        self.z_range = np.arange(zc-dz/2, zc+dz/2, -dens)
        
    def update_atm(self, atm, t=0):
        x, y, z, wx, wy, wz = atm.sample(self.x_range, self.y_range, self.z_range, t)
        nx, ny, nz = x.shape
        points = np.zeros((nx*ny*nz, 4), dtype=np.float32)
        _idx=0
        for ix in range(nx):
             for iy in range(ny):
                 for iz in range(nz):
                     _ival = wz[ix, iy, iz] if np.abs(wz[ix, iy, iz]) >= 0.1 else np.float('nan')
                     points[_idx] = x[ix, iy, iz], y[ix, iy, iz], z[ix, iy, iz] , _ival
                     _idx+=1
        mask = np.isfinite(points[:,3])
        points = points[mask]
        
        self.pc2_msg = sensor_msgs.point_cloud2.create_cloud(self.header, self.msg_fields, points)
        
        
    def publish(self, model):
        stamp = rospy.Time.now()
        self.header.stamp = stamp
        SimplePublisher.publish(self, self.pc2_msg)
        # self.msg.data = np.asarray(points, np.float32).tostring()
        # SimplePublisher.publish(self, self.msg)
        

class GuidanceStatusPublisher(SimplePublisher):
    def __init__(self, topic='guidance/status'):
        SimplePublisher.__init__(self, topic, ros_pat.msg.guidance_status, 'simulator') 

    def publish(self, model):
        msg = ros_pat.msg.guidance_status()
        msg.h = -model.sim.fdm.X[p3_fr.SixDOFAeroEuler.sv_z]
        if len(model.sim.Xees) > 0:
            msg.hdot = -model.sim.Xees[-1][p3_fr.SixDOFEuclidianEuler.sv_zd]
        msg.air_vel = model.sim.fdm.X[p3_fr.SixDOFAeroEuler.sv_va]
        msg.phi = np.rad2deg(model.sim.fdm.X[p3_fr.SixDOFAeroEuler.sv_phi])
        msg.theta = np.rad2deg(model.sim.fdm.X[p3_fr.SixDOFAeroEuler.sv_theta])
        msg.air_vel_sp = model.sim.ctl.get_va_sp()
        msg.phi_sp = np.rad2deg(model.sim.ctl.phi_sp())
        msg.theta_sp = np.rad2deg(model.sim.ctl.theta_sp())
        SimplePublisher.publish(self, msg)

    
class OdomPublisher(SimplePublisher):
    def __init__(self, topic='ros_pat/odom', what='N/A', frame='w_ned'):
        SimplePublisher.__init__(self, topic, nav_msgs.msg.Odometry, what) 
        self.frame = frame
        
    def publish(self, model, _stamp):
        msg = nav_msgs.msg.Odometry()
        msg.header.stamp = _stamp
        msg.header.frame_id = self.frame
        x, y, z = p3_u.get_trans(model)
        #_quat = tf.transformations.quaternion_from_euler(0, 0, 0)
        _quat = tf.transformations.quaternion_from_matrix(model)
        msg.pose.pose = geometry_msgs.msg.Pose(geometry_msgs.msg.Point(x, y, z), geometry_msgs.msg.Quaternion(*_quat))
        # set the velocity
        #odom.child_frame_id = "base_link"
        #odom.twist.twist = Twist(Vector3(vx, vy, 0), Vector3(0, 0, vth))
        SimplePublisher.publish(self, msg)

#
# Joint state
#
class JointStatePublisher(SimplePublisher):
    def __init__(self, topic='/ds_glider/joint_states', what='N/A'):
        SimplePublisher.__init__(self, topic, sensor_msgs.msg.JointState, what) 
        self.msg = sensor_msgs.msg.JointState()
        self.msg.name = ['aileron_left_joint', 'aileron_right_joint', 'elevator_left_joint', 'elevator_right_joint', 'flap_left_joint', 'flap_right_joint']
        self.msg.position = np.zeros(len(self.msg.name))
        self.msg.velocity = []
        self.msg.effort = []

    def publish(self, model, _stamp):
        self.msg.header.stamp = _stamp
        self.msg.position = model
        SimplePublisher.publish(self, self.msg)
#
# Pose
#
class PosePublisher(SimplePublisher):
    def __init__(self, topic='/ds_glider/inertial_vel', what='N/A'):
        SimplePublisher.__init__(self, topic, geometry_msgs.msg.Pose, what)
        self.msg = geometry_msgs.msg.Pose()
        
    def publish(self, model):
        SimplePublisher.publish(self, self.msg)

#
# Text Marker
#              
class TextMarkerPublisher(SimplePublisher):
    def __init__(self, topic, ref_frame, scale=1., argb=(1., 0., 1., 0.), what='N/A'):
        SimplePublisher.__init__(self, topic, visualization_msgs.msg.Marker, what)
        #self.marker_pub = rospy.Publisher(topic, visualization_msgs.msg.Marker, queue_size=1)
        self.marker_msg = visualization_msgs.msg.Marker()
        self.marker_msg.header.frame_id=ref_frame
        self.marker_msg.type = visualization_msgs.msg.Marker.TEXT_VIEW_FACING
        p = self.marker_msg.pose.position; p.x, p.y, p.z = 0, 0, -0.2
        o = self.marker_msg.pose.orientation; o.x, o.y, o.z, o.w = 0, 0, 0, 1
        self.marker_msg.scale.z = scale
        c = self.marker_msg.color; c.a, c.r, c.g, c.b = argb
        self.marker_msg.text = 'foo'

    def publish(self, model):
        self.marker_msg.header.stamp = rospy.Time.now()
        self.marker_msg.text = model
        #p = self.marker_msg.pose.position; p.x, p.y, p.z = carrot_pos[0], carrot_pos[1],carrot_pos[2]
        SimplePublisher.publish(self, self.marker_msg)
 

        
#
# periodic node
#
class PeriodicNode:

    def __init__(self, name):
        rospy.init_node(name)
    
    def run(self, freq):
        rate = rospy.Rate(freq)
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass
