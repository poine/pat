# PAT
Python Aerospace Toolbox

Some stuff related to aerospace simulation


[Documentation index](https://poine.github.io/pat)



### Installing

 * clone sources:
 ```
 cd
 git clone https://github.com/poine/pat.git
 ```

Done... (unless you need the [ROS part](https://poine.github.io/pat/ros_intro.html), for example 3D display in rviz)

Make sure your python is able to find the cloned directory. For example on Linux you could use:

```console
poine@nina:~$ export PYTHONPATH=$PYTHONPATH:/home/poine/pat/src
```
You can test that it worked with

```console
poine@nina:~$ python -c "import pat3; print(pat3.__version__)"
```

 * run examples:

 ** fixed wing
 
```console
poine@nina:~/pat$ ./src/pat3/test/fixed_wing/test_02_att_ctl.py
```




