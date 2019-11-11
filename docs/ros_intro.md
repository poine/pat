---
title: PAT ROS
layout: default
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

$$
\newcommand{\vect}[1]{\underline{#1}}                      % vector
\newcommand{\mat}[1]{\mathbf{#1}}                          % matrices
\newcommand{\est}[1]{\hat{#1}}                             % estimate
\newcommand{\err}[1]{\tilde{#1}}                           % error
\newcommand{\pd}[2]{\frac{\partial{#1}}{\partial{#2}}}     % partial derivatives
\newcommand{\transp}[1]{#1^{T}}                            % transpose
\newcommand{\inv}[1]{#1^{-1}}                              % invert
\newcommand{\norm}[1]{|{#1}|}                              % norm
\newcommand{\esp}[1]{\mathbb{E}\left[{#1}\right]}          % expectation
\newcommand{\identity}[0]{\mathbb{I}}                      % identity
$$


## PAT ROS

PAT comes with some ROS aware tools.

### Installing
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

  * Fixed Wing
```console
poine@nina:~$ roslaunch ros_pat sim_glider.launch
```
 * Multirotor
```console
poine@nina:~$ roslaunch ros_pat sim_multirotor.launch
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

