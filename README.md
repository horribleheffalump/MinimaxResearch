# MinimaxResearch

An unscented Kalman filter (UKF) and conditionally minimax nonlinear filter (CMNF) implementation for 
stochastic dynamic system state estimation.

## UKFilter

Unscented Kalman filter algorithm from [1, Ch. 7.3]. In order to improve the performance, an additional procedure for UKF parameter
tuning is proposed. 

## CMNFilter

Conditionally minimax filter from [2].

## TestEnvironments

The filters are compared on a set of different models:
1. static models: carthesian to polar and sphere coordinates transformation;
2. cubic sensor;
3. -- under construction --


[[1]](http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471369985.html) Kalman Filtering and Neural Networks, Edited by Simon Haykin, 2001, John Wiley & Sons, Inc.

[[2]](http://ieeexplore.ieee.org/document/310035/) Pankov A., Bosov A., Conditionally minimax algorithm for nonlinear system state estimation // IEEE Transactions on Automatic Control, vol. 39(8), 1994, 
pp. 1617-1620. DOI: 10.1109/9.310035.
