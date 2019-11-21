---
title: Fixed Wing FDM
layout: default
---

<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
$$
\newcommand{\vect}[1]{\underline{#1}}                      % vector
\newcommand{\transp}[1]{#1^{T}}                            % transpose
$$

## Flight Dynamic Model

* State Vector
 
   $$
\begin{align}
\vect{X} &= \transp{\begin{pmatrix} x & y & z & v & \alpha & \beta & \phi & \theta & \psi & p & q & r \end{pmatrix}} \\
\vect{X} &= \transp{\begin{pmatrix} \vect{P} & \vect{V} & \vect{E} & \vect{\Omega} \end{pmatrix}}
\end{align}
   $$
  
   where:
    
   * $$ \vect{P} = \transp{\begin{pmatrix} x & y & z \end{pmatrix}} $$ is the position of the center of mass in an euclidian space.
   * $$ \vect{V} = \transp{\begin{pmatrix} v & \alpha & \beta \end{pmatrix}} $$ are the airspeed, angle of attack and sideslip.
   * $$ \vect{E} = \transp{\begin{pmatrix}  \phi & \theta & \psi \end{pmatrix}} $$ is the rotation between world and body frame as euler angles.
   * $$ \vect{\Omega} = \transp{\begin{pmatrix} p & q & r \end{pmatrix}} $$ is the rotational velocity of body frame with respect to world expressed in body frame.




### Kinematics

#### Velocities

$$ 
\dot{v} =  bla

$$

```python
Xdot[sv_v] = np.inner(ivel_body, accel_body)/X[sv_v]
u, v, w = ivel_body
ud, vd, wd = accel_body
Xdot[sv_alpha] = (u*wd - w*ud)/(u**2+w**2)
Xdot[sv_beta] = (X[sv_v]*vd - v*Xdot[sv_v]) / X[sv_v] / math.sqrt(u**2+w**2)
```

#### Attitude
```python
Xdot[sv_phi:sv_psi+1] = p3_alg.euler_derivatives(X_euler, X_rvel_body)

def euler_derivatives(eu, om):

  sph = math.sin(eu[e_phi])
  cph = math.cos(eu[e_phi])
  tth = math.tan(eu[e_theta])
  cth = math.cos(eu[e_theta])

  return [ om[r_p] + sph*tth*om[r_q] + cph*tth*om[r_r],
                         cph*om[r_q] -     sph*om[r_r],
                     sph/cth*om[r_q] + cph/cth*om[r_r] ]
```

### Dynamics
```python
 accel_body = 1./P.m*forces_body - np.cross(X_rvel_body, ivel_body)
 ```


```python
 raccel_body = np.dot(P.invI, m_aero_body + m_eng_body - np.cross(X_rvel_body, np.dot(P.I, X_rvel_body)))
 ```

```console
poine@nina:~/pat$ ./src/pat3/test/fixed_wing/test_02_att_ctl.py
```
