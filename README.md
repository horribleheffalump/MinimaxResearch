# MinimaxResearch

An environment for testing various estimation methods in stochastic dynamic systems:

x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, 

y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t

Filters implemented so far:

## UKFilter

Unscented Kalman filter algorithm from [1, Ch. 7.3]. In order to improve the performance, an additional procedure for UKF parameter
tuning is proposed. [Details](https://github.com/horribleheffalump/MinimaxResearch/blob/master/CMNFvsUT/Docs/UKF.GeneratedXmlDoc.md) 

## EKFilter

Extended Kalman filter algorithm.

## CMNFilter

Conditionally minimax nonlinear filter from [2],[3] and its versions. The whole environment is designed for research in the field of CMN filters, hence the name of the repository. UK end EK filters are added for comparative performance evaluation. 

## TestEnvironments

TestEnvironments project contains a set of tools for modelling studies. It allows to generate sample paths, apply filters, calculate estimation error statistics. Parallelization is supported where it is possible. 

This project contains a set of predefined test environments, implementing (among others) the following models:

Static models:

1. cartesian to polar and sphere coordinates transformation;

Discrete dynamic models:

2. cubic sensor;

3. logistic regression;

4. sampled regression;

5. regression with switching observation.

A comprehensive example of test environment usage is provided in CMNFvsUTFTest project.

TargetTrackingTest project implements a maneuvring target tracking model. It also demonstrates how a continuous dynamic model may be studied in the discrete-time test environment.

## References

[[1]](http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471369985.html) Kalman Filtering and Neural Networks, Edited by Simon Haykin, 2001, John Wiley & Sons, Inc.

[[2]](http://ieeexplore.ieee.org/document/310035/) Pankov A., Bosov A., Conditionally minimax algorithm for nonlinear system state estimation // IEEE Transactions on Automatic Control, vol. 39(8), 1994, 
pp. 1617-1620. DOI: 10.1109/9.310035.

[[3]](https://link.springer.com/article/10.1134%2FS0005117918010010) A. V. Borisov, A. V. Bosov, A. I. Kibzun, G. B. Miller, K. V. Semenikhin, The Conditionally Minimax Nonlinear Filtering Method and Modern Approaches to State Estimation in Nonlinear Stochastic Systems // Automation and Remote Control, Vol. 79(1), 2018, pp. 1–11. DOI: 10.1134/S0005117918010010.
