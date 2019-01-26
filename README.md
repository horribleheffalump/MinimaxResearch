# MinimaxResearch

An environment for testing various estimation methods in stichastic dynamic systems.
Filters implemented so far:

## UKFilter

Unscented Kalman filter algorithm from [1, Ch. 7.3]. In order to improve the performance, an additional procedure for UKF parameter
tuning is proposed. [Details](https://github.com/horribleheffalump/MinimaxResearch/blob/master/CMNFvsUT/Docs/UKF.GeneratedXmlDoc.md) 

## EKFilter

Extended Kalman filter algorithm.

## CMNFilter

Conditionally minimax filter from [2].

## TestEnvironments

The solution contains a set of predefined test environments, implementing (among others) the following models:
Static models:
1. cartesian to polar and sphere coordinates transformation;
Discrete dynamic models:
2. cubic sensor;
3. logistic regression;
4. sampled regression;
5. regression with switching observation;
Continuous dynamic models:
6. maneuvring target tracking (in separate project TargetTrackingTest)


[[1]](http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471369985.html) Kalman Filtering and Neural Networks, Edited by Simon Haykin, 2001, John Wiley & Sons, Inc.

[[2]](http://ieeexplore.ieee.org/document/310035/) Pankov A., Bosov A., Conditionally minimax algorithm for nonlinear system state estimation // IEEE Transactions on Automatic Control, vol. 39(8), 1994, 
pp. 1617-1620. DOI: 10.1109/9.310035.
