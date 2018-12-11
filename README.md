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
 ~/pat/src/real_time_sim_guidance.py
```
  * draw the vehicle
```
rviz -d ~/pat/ros_pat/rviz/multirotor.rviz
```

### Run a batch simulation
* 
```
./src/test_04_sim_guidance.py
```

### View a trajectory

