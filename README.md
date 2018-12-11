# PAT
Python Aerospace Toolbox

A big name for what is now only a simple multirotor flight dynamic model


[Documentation index](https://poine.github.io/pat)



### Installing

 * clone sources:
 ```
 git clone https://github.com/poine/ros_pat.git
 ```
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


### Run a real time simulation
 * run the simulation
```
 ~/pat/src/real_time_sim_guidance.py
```
  * draw the vehicle
```
rviz -d rviz -d ~/pat/ros_pat/rviz/multirotor.rviz
```

### Run a batch simulation
* 
```
./src/test_04_sim_guidance.py
```
