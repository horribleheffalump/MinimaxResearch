# TestEnvironments #

## Type BasicFilter

 Basic class for filters' standardization Siutable for filters with recursive structure, where the estimate on the current step depend only on the current observations, estimate on the previous step, and estimated error covariance 



---
#### Method BasicFilter.Initialize

 Abstract method for filter initialization 



---
#### Method BasicFilter.InitializeAndTrain

 Abstract method for filter initialization and parameter calculation 



---
#### Method BasicFilter.Step(System.Int32,MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Vector{System.Double},MathNet.Numerics.LinearAlgebra.Matrix{System.Double})

 Abstract method for estimate calculation 

|Name | Description |
|-----|------|
|t: |Current time instant|
|y: |Current observation vector|
|xHat: |Estimate on the previous step|
|kHat: |Error covariance matrix|
**Returns**: 



---
#### Method BasicFilter.SaveParamsText

 Virtual method to save parameters as text 



---
#### Method BasicFilter.SaveParams

 Virtual method to save filter parameters as binary for future use 



---
#### Method BasicFilter.LoadParams

 Virtual method to load filter parameters 



---
## Type Properties.Resources

 A strongly-typed resource class, for looking up localized strings, etc. 



---
#### Property Properties.Resources.ResourceManager

 Returns the cached ResourceManager instance used by this class. 



---
#### Property Properties.Resources.Culture

 Overrides the current thread's CurrentUICulture property for all resource lookups using this strongly typed resource class. 



---
#### Property Properties.Resources.LatexPictureTemplte

 Looks up a localized string similar to \centerline{\includegraphics[width=0.9\textwidth]{%file%}} \caption{%caption%}. 



---
#### Property Properties.Resources.OutputFileNameBinTemplate

 Looks up a localized string similar to {name}_{rnd}.pinfo. 



---
#### Property Properties.Resources.OutputFileNameTemplate

 Looks up a localized string similar to {name}_{type}_{0}.txt. 



---
#### Property Properties.Resources.OutputPictureNameTemplate

 Looks up a localized string similar to {name}_{script}_{0}.pdf. 



---
#### Property Properties.Resources.OutputTypeBulk

 Looks up a localized string similar to bulk. 



---
#### Property Properties.Resources.OutputTypeMany

 Looks up a localized string similar to average. 



---
#### Property Properties.Resources.OutputTypeOneObs

 Looks up a localized string similar to sample_obs. 



---
#### Property Properties.Resources.OutputTypeOneState

 Looks up a localized string similar to sample. 



---
## Type TestEnvironmentVector

 Test environment for filters comparison on discrete stochastic dynamic systems 



---
#### Method TestEnvironmentVector.Initialize(System.Int32,System.Int32,System.String,System.Collections.Generic.List{System.ValueTuple{TestEnvironments.FilterType,System.String}},System.Boolean,System.Boolean)

 Initializes the test environment 

|Name | Description |
|-----|------|
|t: |time interval right margin (number of steps)|
|n: |test set size for filters parameter calculation|
|workingFolder: |working folder to store all the output and to load params from|
|filters: |list of filter types |
|save: |if true, save the calculated filter params|
|load: |if true, load the filter params calculated earlier|


---
#### Method TestEnvironmentVector.GenerateBundleSamples(System.Int32,System.Int32,System.String)

 Generates a bundle of sample paths and saves the state dynamics to files (one file for each coordinate) 

|Name | Description |
|-----|------|
|t: |time interval right margin (number of steps)|
|n: |number of sample paths to generate|
|outputFolder: |folder to save the results|


---
#### Method TestEnvironmentVector.GenerateOne(System.String,System.Nullable{System.Int32})

 Generates a single sample path (bundle of sample paths, if n>1) and saves the state dynamics and observations to separate files (one file for each coordinate of each sample) 

|Name | Description |
|-----|------|
|outputFolder: |folder to save the results|
|n: |number of sample paths to generate|


---
#### Method TestEnvironmentVector.GenerateBundles(System.Int32,System.Int32,System.String,System.Boolean,System.Int32,System.Boolean,System.Boolean)

 Generates N bundles of sample paths and applies the filters. The statistics for estimate errors is calculated by taking average on each bundle, and then on the whole set of bundles. 

|Name | Description |
|-----|------|
|N: |Number of bundles|
|n: |Number of samples in each bundle|
|outputFolder: |folder to save the results|
|parallel: |if true, bundles are calculated in parallel|
|parallel_degree: |degree of parallelizm|
|doSaveBin: |if true, the results are saved in binary format (allows subsequent aggregation)|
|doSaveText: |if true, the results are saved in text format (subsequent aggregation is not possible)|


---
#### Method TestEnvironmentVector.CalculateBundle(System.Int32,TestEnvironments.ProcessInfo)

 Generates a bundle of sample paths, applies the filters and returns estimate errors statistics 

|Name | Description |
|-----|------|
|n: |Number of samples in a bundle|
|processInfo: |Structure with initial parameters, e.g. filter list|
**Returns**: Estimate error statistics for generated bundle



---
#### Method TestEnvironmentVector.Aggregate(System.String,System.String,System.Boolean,System.Boolean)

 Imports and aggregates data from prevoiously generated bundles of trajectories. The statistics for estimate errors is calculated by taking average the whole set. 

|Name | Description |
|-----|------|
|inputFolder: |folder to search files with previously generated statistics subject to aggregation|
|outputFolder: |folder to save the results|
|doSaveBin: |if true, the results are saved in binary format (allows subsequent aggregation)|
|doSaveText: |if true, the results are saved in text format (subsequent aggregation is not possible)|


---
#### Method TestEnvironmentVector.RunScript(System.String,System.String)

 Runs a Python script with single param 

|Name | Description |
|-----|------|
|scriptName: |script file path|
|scriptParam: |script parameter as string|


---


