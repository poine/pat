#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


x = np.array([0, 1, 0])
y = np.array([0, 1, 2])

x = np.array([ 0. ,  1.2,  1.9,  3.2,  4. ,  6.5])
y = np.array([ 0. ,  2.3,  3. ,  4.3,  2.9,  3.1])

t, c, k = interpolate.splrep(x, y, s=0, k=4)
print(f'''\
t: {t}
c: {c}
k: {k}
''')
N = 100
xmin, xmax = x.min(), x.max()
xx = np.linspace(xmin, xmax, N)
spline = interpolate.BSpline(t, c, k, extrapolate=False)

plt.plot(x, y, 'bo', label='Original points')
plt.plot(xx, spline(xx), 'r', label='BSpline')

plt.show()
