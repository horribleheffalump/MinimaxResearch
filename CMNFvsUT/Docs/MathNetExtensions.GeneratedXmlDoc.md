# MathNetExtensions #

#### Method Exts.Inverse(MathNet.Numerics.LinearAlgebra.Matrix{System.Double},System.Double,System.Double,System.Int32)

 Inverts a matrix using Jacobi method 

|Name | Description |
|-----|------|
|x: |Square matrix to invert|
|tol1: |upper bound of off - diagonal elements (default 1e-32)|
|tol2: |lower bound of diagonal non-zero elements (default 1e-32)|
**Returns**: 



---
## Type FiniteDiscreteDistribution

 Random variable generator for a finite discrete distribution given by its probability measure 



---
## Type RandomOptimizer

 Random shoot optimizer: inputs of the objective function are randomly sampled and the best sample is chosen as optimal 



---
#### Method RandomOptimizer.Minimize(System.Func{MathNet.Numerics.LinearAlgebra.Vector{System.Double},System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},System.Int32,System.Int32,System.String)

Random shoot optimization procedure.

- Step 1: generate PointsUniform uniform random samples in [LowerBound, UpperBound] interval, calculate the objective outputs.

- Step 2: choose the random sample with the best output and generate PointsNormal random uniform samples in closer area of this best sample.

- Step 3: again choose the random sample with the best output, save the samples ordered by the objective value to the output file and return the best found sample and objective output.

Samples generation and objective function calculation is performed asynchronously.



|Name | Description |
|-----|------|
|Objective: |Objective function|
|LowerBound: |Lower bound of the uniform sampling interval|
|UpperBound: |Upper bound of the uniform sampling interval|
|PointsUniform: |Number of uniform samples on the first step|
|PointsNormal: |Number of normal samples on the second step (may be zero)|
|OutputFileName: |Output file name|
**Returns**: Returns touple (best found objective value, sample with the best found objective value)



---
## Type RandomVector`1

 Random vector with given expextation and covariance matrix generator. 

|Name | Description |
|-----|------|
|T: |Continuous distribution|


---


