# CMNF #

## Type BCMNFilter

A version of Conditionnaly Minimax Nonlinear filter for a model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)

Usage:

initialize with EstimateParameters method, which calculates the filter coefficients using the provided test set.

After the filter is initialized, it may be used for estimate step-by-step calculation

This version takes into account some additional correlations





---
## Type ModifiedCMNFilter

Modified version of Conditionnaly Minimax Nonlinear filter for a model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)

In contrast with original CMNF does not need initialization. The filter coefficients are calculated online.





---
## Type CMNFilter

Original Conditionnaly Minimax Nonlinear filter for a model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)

Usage:

initialize with EstimateParameters method, which calculates the filter coefficients using the provided test set.

After the filter is initialized, it may be used for estimate step-by-step calculation





---


