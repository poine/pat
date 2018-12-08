---
title: Multirotor FDM
layout: default
---


<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>


Let us consider a vehicle comprising a set of $$N$$ identical power trains $$R_i, \quad i \in[1:N]$$ located at coordinates $$(X_i,Y_i), i\in[1:N]$$ and spinning in the same plane in the direction $$D_i, i\in[1:N], D_i\in[-1;1]$$ at a rotational speed $$\omega_i, i\in[1:N]$$. 


Assuming a quasi static regime, the force produced by each rotor can be considered normal to the rotor plane and proportional to the square of its rotational speed. Under the same asumption, the torque produced by each rotor can also be assumed to be in the same direction and proportional to the square of the rotational speed.

$$
\begin{equation}
\vect{M}_{i}^{B} = \begin{pmatrix}X_i C_t \omega_i^2\\Y_i C_t \omega_i^2\\D_i C_m \omega_i^2\end{pmatrix}
\end{equation}
$$

Where $$C_t$$ is a **thrust** coefficient and $$Cd$$ is a **torque** coefficient. It has been shown experimentaly that $$\frac{Ct}{Cm} \approx 10$$ on the mikrokopter.
