# EKF #

## Type ExtendedKalmanFilter

Extended Kalman filter for a model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)

Usage:

Provide system dynamics and observations functions Phi1, Phi2, Psi1, and Psi2 as well as their derivatives dPhi and dPsi.

If the derivatives are not provided, the original function is considered linear and its derivative is a constant matrix which is calculated columnwise.

A particular funcion may be provided for the prediction and its covariance calculation. That is made for the continuous dynamics case, where predicted covariance estimate kTilde is calculated as a solution to the Riccati equation





---


