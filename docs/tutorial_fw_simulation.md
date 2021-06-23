---
title: Fixed Wing Simulation
layout: default
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

$$
X = \Omega
$$

# Fixed Wing Aircraft Simulation

[Simulating](https://en.wikipedia.org/wiki/Computer_simulation) an aircraft means computing its evolutions through time. It can be achieved by representing the aircraft as a solid, computing the [forces and moments](https://en.wikipedia.org/wiki/Aerodynamics) applied to it, and (numerically) solving the differential equations resulting for the [solid's physics](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion). Pat provides tools to easily simulate an aircraft using the Python language.



## Quick start:

```
import numpy as np, pat3.utils as p3_u, pat3.vehicles.fixed_wing.legacy_6dof as p3_fw
dm = p3_fw.DynamicModel_ee(p3_u.pat_ressource('data/vehicles/cularis.xml'))
Xe, Ue = dm.trim({'h':0, 'va':10, 'gamma':0}, report=True)
time = np.arange(0, 5., 0.01)
X, U = np.zeros((len(time), dm.sv_size)), np.zeros((len(time), dm.iv_size))
X[0] = dm.reset(X0, time[0])
for i in range(1, len(time)):
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1])
dm.plot_trajectory(time, X, U)
plt.show()
```



```
./src/pat3/test/fixed_wing/test_01_dynamics.py
```

