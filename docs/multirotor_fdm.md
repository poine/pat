---
title: Multirotor FDM
layout: default
---

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
   * $$ \vect{\Omega} = \transp{\begin{pmatrix} p & q & r \end{pmatrix}} $$ is the rotationale velocity of body frame with respect to world expressed in body frame.


 * Input Vector

   $$
\vect{U} = \transp{\begin{pmatrix} \vect{F}_b & \vect{M}_b \end{pmatrix}}
   $$

   where:
   
    * $$ \vect{F}_b $$ is the sum of external forces applied to the solid, expressed in body frame
    * $$ \vect{M}_b $$ is the total moment created by external forces,  expressed in body frame

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

```py
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


Let us consider a vehicle comprising a set of $$N$$ identical power trains $$R_i,  i \in[1:N]$$ located at coordinates $$(X_i,Y_i), i\in[1:N]$$ and spinning in the same plane in the direction $$D_i, i\in[1:N], D_i\in[-1;1]$$ at a rotational speed $$\omega_i, i\in[1:N]$$. 


Assuming a quasi static regime, the force produced by each rotor can be considered normal to the rotor plane and proportional to the square of its rotational speed. Under the same asumption, the torque produced by each rotor can also be assumed to be in the same direction and proportional to the square of the rotational speed.

$$
\begin{equation}
\vect{M}_{i}^{B} = \begin{pmatrix}X_i C_t \omega_i^2\\Y_i C_t \omega_i^2\\D_i C_m \omega_i^2\end{pmatrix}
\end{equation}
$$

Where $$C_t$$ is a **thrust** coefficient and $$Cd$$ is a **torque** coefficient. It has been shown experimentaly that $$\frac{Ct}{Cm} \approx 10$$ on the mikrokopter.
