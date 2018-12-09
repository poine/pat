---
title: Multirotor FDM
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


## Solid Dynamics

 * State Vector
 
   $$
\begin{align}
\vect{X} &= \transp{\begin{pmatrix} x & y & z & \dot{x} & \dot{y} & \dot{z} & qi & qx & qy & qz & p & q & r \end{pmatrix}} \\
\vect{X} &= \transp{\begin{pmatrix} \vect{P} & \vect{V} & \vect{Q} & \vect{\Omega} \end{pmatrix}}
\end{align}
   $$
  
   where:
    
   * $$ \vect{P} = \transp{\begin{pmatrix} x & y & z \end{pmatrix}} $$ is the position of the center of mass in an euclidian space.
   * $$ \vect{V} = \transp{\begin{pmatrix} \dot{x} & \dot{y} & \dot{z} \end{pmatrix}} $$ is the velocity os the CG in world frame.
   * $$ \vect{Q} = \transp{\begin{pmatrix} qi & qx & qy & qz \end{pmatrix}} $$ is the rotation between world and body frame as a unit quaternion.
   * $$ \vect{\Omega} = \transp{\begin{pmatrix} p & q & r \end{pmatrix}} $$ is the rotational velocity of body frame with respect to world expressed in body frame.


 * Input Vector

   $$
\vect{U} = \transp{\begin{pmatrix} \vect{F}_b & \vect{M}_b \end{pmatrix}}
   $$

   where:
   
    * $$ \vect{F}_b $$ is the sum of external forces applied to the solid, expressed in body frame.
    * $$ \vect{M}_b $$ is the total moment created by external forces,  expressed in body frame.

### Kinematics

  *  $$\dot{\vect{P}}$$ is simply contained in $$\vect{X}$$
  
     $$
\dot{\vect{P}} = \vect{V}
	 $$ 
  *  $$\dot{\vect{Q}}$$ is computed from  $$\vect{Q}$$ and $$\vect{\Omega}$$ as:
  
     $$ 
\dot{\vect{Q}} = -\frac{1}{2} \text{skewsym}(\vect{\Omega}) \vect{Q}
	 $$

### Translational Dynamics

$$
\dot{\vect{V}} = \frac{1}{m} \left( F^w + \transp{\begin{pmatrix} 0 & 0 & mg \end{pmatrix}} \right)
$$



### Rotational Dynamics

$$
\dot{\vect{\Omega}} = \inv{\mat{J}} \left( \vect{M}^b - \vect{\Omega} \wedge \mat{J} \vect{\Omega} \right)
$$


### Code

```python
def solid_cont_dyn(X, F_b, M_b, P):
    Xd = np.zeros(sv_size)
    p_w, v_w, q_w2b, om_b = X[sv_slice_pos], X[sv_slice_vel], X[sv_slice_quat], X[sv_slice_rvel]
    # Translational kinematics
    Xd[sv_slice_pos] = v_w
    # Newton for forces
    R_w2b =  pal.rmat_of_quat(q_w2b)
    Xd[sv_slice_vel] = 1./P.m*(np.dot(R_w2b.T, F_b) + [0, 0, P.m*P.g])
    # Rotational kinematics
    Xd[sv_slice_quat] = pal.quat_derivative(q_w2b, om_b)
    # Newton for moments
    Xd[sv_slice_rvel] = np.dot(P.invJ, M_b - np.cross(om_b, np.dot(P.J, om_b)))
    return Xd
```

[link](https://github.com/poine/pat/blob/master/src/pat3/vehicles/rotorcraft/multirotor_fdm.py)


## Multirotor

