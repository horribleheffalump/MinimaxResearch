# UKF #

## Type SigmaPoints

 Static sigma points generator 



---
#### Method SigmaPoints.Generate(MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double)

Generates sigma-points given the mean x, covarianve P and spread parameter lambda

- Xi_0 = x

- Xi_i = x + sqrt((L+lambda)P)_i, i = 1,...,L

- Xi_i = x - sqrt((L+lambda)P)_{i-L}, i = L+1,...,2L 

L - dimention of x, sqrt(P) - Cholesky decomposition 



|Name | Description |
|-----|------|
|x: |Mean|
|P: |Covariance|
|lambda: |Spread parameter|
**Returns**: 



---
## Type UnscentedTransform

 The unscented transform 



---
#### Method UnscentedTransform.Transform(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},UKF.UTParams,MathNet.Numerics.LinearAlgebra.Vector{System.Double}@,MathNet.Numerics.LinearAlgebra.Matrix{System.Double}@,MathNet.Numerics.LinearAlgebra.Matrix{System.Double}@)

 The unscented transform 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|mX: |Mean of the transformed random variable|
|dX: |Cov of the ransformed random variable|
|dY: |Cov of the additive random variable|
|p: |Parameters of the unscented transform|
|y: |Returns: approximated mean of the transformed variable|
|Kxy: |Returns: approximated cross-covariance of the initial and the transformed variable|
|Kyy: |Returns: approximated covariance of the transormed variable|


---
## Type UTDefinitionType

 Unscented transform parameters optimization type. 



---
## Type UTParams

Parameters of the unscented transform:

- Lambda - scaling parameter

- Wm - weights for sample mean

- Wc - weights for sample covariance

Three ways to define:

- explicitly,

- implicitly with one parameter alpha0,

- implicitly with three parameters alpha, beta, kappa





---
#### Method UTParams.#ctor(System.Int32,System.Double[])

 Constructor 

|Name | Description |
|-----|------|
|L: |Dimension of the transformed random variable|
|p: |Input params to define the unscented transform params. The appropriate method to define UT params is chosen depending on the size of the array.|


---
#### Method UTParams.SetUTParams(System.Int32,System.Double)

 The parameters of the unscented transform are defined by a single parameter alpha0: 

|Name | Description |
|-----|------|
|L: |Dimension of the transformed random variable|
|alpha0: |Alpha0 - weight of the central points for both sample mean and cov |


---
#### Method UTParams.SetUTParams(System.Int32,System.Double,System.Double,System.Double)

 The parameters of the unscented transform are defined by three parameters: alpha, beta, kappa 

|Name | Description |
|-----|------|
|L: |Dimension of the transformed random variable|
|alpha: |Alpha - determines the spread of the sigma points around the mean of the transformed random variable (small positive value 0 \leq alpha \leq 10e-4)|
|beta: |Beta - is used to incorporate prior knowledge of the distribution of the transformed random variable (for Gaussian b = 2 is optimal)|
|kappa: |Kappa - is a secondary scaling parameter|


---
#### Method UTParams.SetUTParams(System.Int32,System.Double,System.Double,System.Double,System.Double)

 Explicit definition of the unscented transform parameters 

|Name | Description |
|-----|------|
|L: |Dimension of the transformed random variable|
|lambda: |Scaling parameter|
|wm0: |Central point weight for the sample mean|
|wc0: |Central point weight for the sample cov|
|wi: |Non-central points weight for sample mean and cov|


---
## Type UTStaticEstimate

Unscented transform estimate for the model y = Phi(x) + Nu, where

- x is a random variable with known mean and covariance,

- Nu - noise with zero mean and known covariance

Usage:

- specify the parameters of the unscented transform in utProperty manually or by means of primitive optmization procedure: UTStaticEstimate.EstimateParameters

- calculate the estimate: UTStaticEstimate.Estimate

It should be noted, that the train and test sets may be different. That is, the arrays of samples X = [x_0,...,x_N] and observations Y = [y_0,...,y_N] = [Phi(x_0) + Nu_0,...,Phi(x_N) + Nu_N] may vary for the step of the unscented transform parameters optimization and the step of unscented transform estimate calculation.





---
#### Method UTStaticEstimate.#ctor(UKF.UTDefinitionType)

 Constructor 

|Name | Description |
|-----|------|
|type: |The way to calculate the unscented transform parameters in the optimization procedure.|


---
#### Method UTStaticEstimate.Estimate(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}@,MathNet.Numerics.LinearAlgebra.Matrix{System.Double}@,MathNet.Numerics.LinearAlgebra.Matrix{System.Double}@)

 Provides the unscented transform estimate for the array X = [x_0,...,x_N] of random variable x samples given the array of observations Y = [y_0,...,y_N], where y_i = Phi(x_i) + Nu_i. Parameters of the unscented transform utParams should be initialized. 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
|mErr_UT: |Returns: estimation error mean vector|
|KErr_UT: |Returns: estimation error covariance marix|
|KErrTh_UT: |Returns: estimation error theoretical covariance marix|
**Returns**: Array of estimates \hat{X} = [\hat{x}_0,...,\hat{x}_N]



---
#### Method UTStaticEstimate.EstimateParameters(System.Int32,System.Int32,System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.String)

 Calls the static unscented transform parameters "optimization" procedure and saves the result into the utParams property. The way how the random samples are transformed to the UT params is determined by the optimizationType property. 

|Name | Description |
|-----|------|
|N1: |Number of samples on the step 1|
|N2: |Number of samples on the step 2|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|Crit: |Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat)) |
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
|outputFolder: |The results are saved to this folder in file "UT_optimization_{type}.txt"|


---
#### Method UTStaticEstimate.UTParmsOptimize(System.Int32,System.Int32,UKF.UTDefinitionType,System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.String)

Unscented transform parameters "optimization" procedure.

- Step 1: generate N1 random samples, calculate the unscented transform estimates given the parameters determined by each random sample.

- Step 2:choose the random sample with the best estimate quality criterion value and generate N2 random samples in closer area of this sample.

- Step 3:again choose the random sample with the best estimate quality criterion value, save the samples ordered by the criterion value to the output file and return the best found unscented transform parameters.

The UTOptimizationType type param determines the way how the random samples define the unscented tranform params (UTParams).

- If type is UTOptimizationType.ImplicitAlpha, then random samples define alpha0 - scalar weight of the central points for both sample mean and cov: UTParams.SetUTParams(int, double);

- If type is UTOptimizationType.ImplicitAlphaBetaKappa, then random samples are vectors of dim 3 and represent three parameters alpha, beta, kappa which are then transformed to the parameters of the inscented transform: UTParams.SetUTParams(int, double, double, double);

- If type is UTOptimizationType.Explicit, then random samples are vectors of dim 4 and explicitly define the unscented transform parameters: UTParams.SetUTParams(int, double, double, double, double). ///TODO it is not right to define the parameters of the unsctnted transform arbitraty, they have to be interdependent, so that the mean and cov would be transformed correctly.



|Name | Description |
|-----|------|
|N1: |Number of samples on the step 1|
|N2: |Number of samples on the step 2|
|type: |The way how the random samples are transformed to the UT params|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|Crit: |Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat)) |
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
|outputFolder: |The results are saved to this folder in file "UT_optimization_{type}.txt"|
**Returns**: Returns the parameters of the unscented transform with best estimation quality



---
#### Method UTStaticEstimate.CalculateSampleCriterion(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.Distributions.IContinuousDistribution[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Generates a random sample for the unscented transform parameters and calculates the criterion value for the unscented transform estimate. The size of the distribution parameter determines the unscented transform parameters definition method (UTParams) 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|Crit: |Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat)) |
|distribution: |Array of distributions to generate random unscented transform parameters|
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
**Returns**: 



---
#### Method UTStaticEstimate.CalculateCriterionValue(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},UKF.UTParams,MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Calculates the criterion value for the estimate given the particular unscented transform parameters 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|Crit: |Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat)) |
|p: |Unscented transform parameters|
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
**Returns**: 



---


