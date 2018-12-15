# PAT
Python Aerospace Toolbox

A big name for what is now only a simple multirotor flight dynamic model


[Documentation index](https://poine.github.io/pat)



### Installing

 * clone sources:
 ```
 cd
 git clone https://github.com/poine/pat.git
 ```

Done... unless you need the ROS part (real time display in rviz), then

* install ROS and create a catkin workspace
 ```
mkdir -p ~/catkin_ws/src
cd ~/catkin_ws/src
catkin_init_workspace
cd ~/catkin_ws
catkin_make
```
  * create symlink to ros package directory in workspace
```
cd ~/catkin_ws/src
ln -s ~/pat/ros_pat .
``` 
  * build 
```
cd ~/catkin_ws/
catkin_make
```
 * and of course source the workspace config
 ```
source ~/catkin_ws/src/devel/setup.bash
```

### Run a real time simulation
 * run the simulation
```
rosrun ros_pat sim_guidance.py _traj_name:=oval_with_intro _time_factor:=0.5
```
  * view the simulation in rviz
```
rviz -d ~/pat/ros_pat/rviz/multirotor.rviz
```

### Run a batch simulation
* 
```
./src/sim_guidance.py
```

### View a trajectory

```
rosrun ros_pat display_trajectory.py _traj_name:=oval_with_intro
```

```
rosrun ros_pat display_track.py
```

```
./src/plot_trajectory.py --traj smooth_back_and_forth
```
