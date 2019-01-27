# UKF #

## Type UKFilter

Unscented Kalman filter for a model x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t, y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t, W_t ~ (M_w, R_w), Nu_t ~ (M_nu, R_nu)

Usage:

- specify the parameters of the unscented transform on forecast and correction phases in utParamsForecast and utParamsCorrection properties manually or by means of the optmization procedure: UKFilter.EstimateParameters

- calculate the estimate step by step with UKFilter.Step

The train and test trajectory bundles may be different. That is, the arrays of discrete vector models may vary for the step of the unscented transform parameters optimization and the step of unscented transform filter calculation.





---
#### Method UKFilter.#ctor(UKF.UTDefinitionType,UKF.OptimizationMethod)

 Constructor 

|Name | Description |
|-----|------|
|type: |Unscented transform parameters definition type|
|method: |Unscented transform parameters optimization method|


---
#### Method UKFilter.EstimateParameters(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},System.Int32,NonlinearSystem.DiscreteVectorModel[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.String)

 Calls the static unscented transform parameters optimization procedure and saves the result into the utParamsForecast and utParamsCorrection properties. The way how to define the UT params is determined by the optimizationType property. The optimization method is determined by the optimizationMethod property. 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|Crit: |Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T)) |
|T: |The upper bound of the observation interval|
|models: |Discrete vector model samples|
|xhat0: |Initial condition|
|DX0Hat: |Initial condition covariance|
|outputFolder: |The results are saved to this folder in file "UT_optimization_{type}.txt"|


---
#### Method UKFilter.EstimateParametersStepwise(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},System.Int32,NonlinearSystem.DiscreteVectorModel[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.String)

 Calls the static unscented transform parameters stepwise optimization procedure and saves the result into the utParamsForecast and utParamsCorrection properties. The way how to define the UT params is determined by the optimizationType property. The optimization method is determined by the optimizationMethod property. 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|Crit: |Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T)) |
|T: |The upper bound of the observation interval|
|models: |Discrete vector model samples|
|xhat0: |Initial condition|
|DX0Hat: |Initial condition covariance|
|outputFolder: |The results are saved to this folder in file "UT_optimization_{type}.txt"|


---
#### Method UKFilter.DefineOptimizationParameters(UKF.UTDefinitionType,MathNet.Numerics.LinearAlgebra.Vector{System.Double},System.String)

 Defines the number of dimensions, lower and upper bound and the initial guess and the file name to store results for the Unscented transform parameters optimization procedures. 

|Name | Description |
|-----|------|
|type: |Unscented transform parameters definition type|
|xhat0: |Initial condition|
|filename: |filename name template to substitute the optimization type name"|
**Returns**: 



---
#### Method UKFilter.SampleVectorToUTParams(MathNet.Numerics.LinearAlgebra.Vector{System.Double},System.Int32)

 The sample vector is transformed to a couple of the unscented transform parameters (the way it is done depends on the length of the sample vector, see UTParams and its constructors for details). 

|Name | Description |
|-----|------|
|P: |Sample vector to be transformed to a couple of the unscented transfrom parameters|
|dim: |The state vector dimention|
**Returns**: 



---
#### Method UKFilter.CalculateSampleCriterion(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},System.Int32,NonlinearSystem.DiscreteVectorModel[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Wrapper for CalculateCriterionValue function. Calculates the criterion value given the provided unscented transform parameters for the forecast and the correction phases. 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|Crit: |Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T)) |
|P: |Sample vector to be transformed to a couple of the unscented transfrom parameters|
|T: |The upper bound of the observation interval|
|models: |Discrete vector model samples|
|xhat0: |Initial condition|
|DX0Hat: |Initial condition covariance|
**Returns**: The criterion value for the unscented transfrom parameters obtained from the sample vactor



---
#### Method UKFilter.CalculateSampleStepwiseCriterion(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},System.Int32,NonlinearSystem.DiscreteVectorModel[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Matrix{System.Double}[])

 Wrapper for CalculateStepwiseCriterionValue function. Calculates the criterion value given the provided unscented transform parameters for the forecast and the correction phases. 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|Crit: |Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T)) |
|P: |Sample vector to be transformed to a couple of the unscented transfrom parameters|
|t: |Step|
|models: |Discrete vector model samples|
|xHat: |array of the state estimates on the previous step for all the given models|
|PHat: |array of the covariance matricies on the previous step for all the given models|
**Returns**: The criterion value for the particular unscented transform parameters



---
#### Method UKFilter.CalculateCriterionValue(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},UKF.UTParams,UKF.UTParams,System.Int32,NonlinearSystem.DiscreteVectorModel[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Calculates the criterion value for the estimate given the particular unscented transform parameters 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|Crit: |Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T)) |
|p1: |Unscented transfrom parameters for the forecast phase|
|p2: |Unscented transfrom parameters for the correction phase|
|T: |The upper bound of the observation interval|
|models: |Discrete vector model samples|
|xhat0: |Initial condition|
|DX0Hat: |Initial condition covariance|
**Returns**: The criterion value for the particular unscented transform parameters



---
#### Method UKFilter.CalculateStepwiseCriterionValue(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},UKF.UTParams,UKF.UTParams,System.Int32,NonlinearSystem.DiscreteVectorModel[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Matrix{System.Double}[])

 Calculates the criterion value for the estimate on a particular step given the particular unscented transform parameters 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|Crit: |Criterion: a function which determines the quality of the unscented Kalman filter. Depends on the sample covariance of the estimation error on the last step: val = Crit(Cov(X_T-Xhat_T,X_T-Xhat_T)) |
|p1: |Unscented transfrom parameters for the forecast phase|
|p2: |Unscented transfrom parameters for the correction phase|
|t: |Step|
|models: |Discrete vector model samples|
|xHat: |array of the state estimates on the previous step for all the given models|
|PHat: |array of the covariance matricies on the previous step for all the given models|
**Returns**: The criterion value for the particular unscented transform parameters



---
#### Method UKFilter.Step(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},UKF.UTParams,UKF.UTParams,System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Performs a step of Unscented Kalman Filter given the particular unscented transform parameters for forecast and correction phases 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|p1: |Unscented transfrom parameters for the forecast phase|
|p2: |Unscented transfrom parameters for the correction phase|
|t: |Current step time instant|
|y: |Observations on the current step|
|xHat_: |Estimate on the previous step|
|P_: |Approximated previous step error covariance|
|xHat: |Returns: current step estimate|
|P: |Returns: approximated current step error covariance|


---
#### Method UKFilter.Step(System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Performs a step of Unscented Kalman Filter with fixed transformation parameters for forecast and correction phases (utParamsForecast and utParamsCorrection must be initialized) 

|Name | Description |
|-----|------|
|Phi1: |State transformation: a nonlinear function which determines the dynamics: x_{t+1} = Phi_1(x_t) + Phi_2(x_t) W_t|
|Phi2: |Noise multiplicator in the dynamics equation: x_{t+1} = Phi(x_t) + W_t|
|Psi1: |Observations transformation: a nonlinear function which determines the relation between the state and the observations: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Psi2: |Noise multiplicator in the observations equation: y_t = Psi_1(x_t) + Psi_2(x_t) Nu_t|
|Mw: |Mean of the noise in the dynamics equation |
|Rw: |Covariance matrix of the state disturbances|
|Mnu: |Mean of the noise in the obseration equation |
|Rnu: |Convariance matrix of the observation noise|
|t: |Current step time instant|
|y: |Observations on the current step|
|xHat_: |Estimate on the previous step|
|P_: |Approximated previous step error covariance|
|xHat: |Returns: current step estimate|
|P: |Returns: approximated current step error covariance|


---
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
**Returns**: Matrix of sigma points



---
## Type UnscentedTransform

 The unscented transform 



---
#### Method UnscentedTransform.Transform(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},UKF.UTParams,MathNet.Numerics.LinearAlgebra.Vector{System.Double}@,MathNet.Numerics.LinearAlgebra.Matrix{System.Double}@,MathNet.Numerics.LinearAlgebra.Matrix{System.Double}@)

 The unscented transform for y = Phi(x) + nu 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|mX: |Mean of the transformed random variable|
|dX: |Cov of the transformed random variable|
|dNu: |Cov of the additive random variable|
|p: |Parameters of the unscented transform|
|y: |Returns: approximated mean of the transformed variable|
|Kxy: |Returns: approximated cross-covariance of the initial and the transformed variable|
|Kyy: |Returns: approximated covariance of the transormed variable|


---
## Type OptimizationMethod

 Unscented transform parameters optimization method. 



---
## Type UTDefinitionType

 Unscented transform parameters definition type. 



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
#### Property UTParams.Params

 Get unscented transformation paramters as array 



---
## Type UTStaticEstimate

Unscented transform estimate for the model y = Phi(x) + Nu, where

- x is a random variable with known mean and covariance,

- Nu - noise with zero mean and known covariance

Usage:

- specify the parameters of the unscented transform in utProperty manually or by means of the optmization procedure UTStaticEstimate.EstimateParameters

- calculate the estimate: UTStaticEstimate.Estimate

It should be noted, that the train and test sets may be different. That is, the arrays of samples X = [x_0,...,x_N] and observations Y = [y_0,...,y_N] = [Phi(x_0) + Nu_0,...,Phi(x_N) + Nu_N] may vary for the step of the unscented transform parameters optimization and the step of unscented transform estimate calculation.





---
#### Method UTStaticEstimate.#ctor(UKF.UTDefinitionType,UKF.OptimizationMethod)

 Constructor 

|Name | Description |
|-----|------|
|type: |Unscented transform parameters definition type|
|method: |Unscented transform parameters optimization method|


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
#### Method UTStaticEstimate.EstimateParameters(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.String)

 Calls the static unscented transform parameters optimization procedure and saves the result into the utParams property. The way how to define the UT params is determined by the optimizationType property. The optimization method is determined by the optimizationMethod property. 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|Crit: |Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat)) |
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
|outputFolder: |The results are saved to this folder in file "UT_optimization_{type}.txt"|


---
#### Method UTStaticEstimate.CalculateSampleCriterion(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}},System.Func{MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double}[],MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Wrapper for CalculateCriterionValue function. The sample vector is transformed to the form unscented transform parameters (the way it is done depends on the length of the sample vector, see UTParams and its constructors for details). Then the criterion value given the provided unscented transform parameters is calculated. 

|Name | Description |
|-----|------|
|Phi: |Transformation: a nonlinear function which determines the transformation of the random vector variable: y = Phi(x) + nu|
|Crit: |Criterion: a function which determines the quality of the unscented transform estimate. Depends on the sample covariance of the estimation error: val = Crit(Cov(X-Xhat,X-Xhat)) |
|P: |Sample vector to be transformed ro unscented transfrom parameters|
|X: |Array of initial variable x samples|
|Y: |Array of transformed variable y = Phi(x) + nu samples|
|mX: |Mean of x|
|KX: |Cov of x|
|KNu: |Cov of the noize nu|
**Returns**: The criterion value for the unscented transfrom parameters obtained from the sample vactor



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
**Returns**: The criterion value for the particular unscented transform parameters



---


