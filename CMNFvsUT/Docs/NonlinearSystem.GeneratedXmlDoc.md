# NonlinearSystem #

## Type DiscreteVectorModel

Model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)

Usage:

Provide starting point X0, system dynamics and observations functions Phi1, Phi2, Psi1, Psi2, and noise generators W and Nu.

Model path sample is generated step-by-step and saved in a dictionary (t, Vector(State,Observations)).





---


