# TestEnvironments #

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

 Test environment for CMN and UT filters comparison on discrete stochastic dynamic systems 



---
#### Method TestEnvironmentVector.Initialize(System.Int32,System.Int32,System.Int32,System.String,System.Collections.Generic.List{TestEnvironments.FilterType})

 Initializes the test environment by calculating the statistics for CMN and UT filters 

|Name | Description |
|-----|------|
|t: |time interval right margin (number of steps)|
|n: |number of trajectories|
|doCalculateUKF: ||
|doCalculateUKFStepwise: ||
|outputFolder: ||


---
#### Method TestEnvironmentVector.GenerateBundleSamples(System.Int32,System.Int32,System.String)

 Generates a bundle of trajectories and saves the state dynamics to files 

|Name | Description |
|-----|------|
|t: |time interval right margin (number of steps)|
|n: |number of trajectories|


---
#### Method TestEnvironmentVector.GenerateBundle(System.Int32,System.String)

 Generates a bundle of trajectories, applies the CMN and UK filters, calculates statistics for estimate errors and saves it to files 

|Name | Description |
|-----|------|
|n: |Number of trajectories|
|folderName: |Output folder name|
|doCalculateUKF: |(optional, default = true) if true, UKF and CMNF estimates are calculated, if false - only CMNF |


---
#### Method TestEnvironmentVector.GenerateBundles(System.Int32,System.Int32,System.String,System.Boolean)

 Generates N bundles of trajectories, applies the CMN and UK filters. The statistics for estimate errors are calculated by taking average on each bundle, and then on the whole set of bundles. Same as GenerateBundle but for larger numbers of trajectories. 

|Name | Description |
|-----|------|
|N: |Number of bundles|
|n: |Number of trajectories|
|folderName: |Output folder name|
|doCalculateUKF: |(optional, default = true) if true, UKF and CMNF estimates are calculated, if false - only CMNF |


---
#### Method TestEnvironmentVector.RunScript(System.String,System.String[])

 Runs a python script to process the calculated test data 

|Name | Description |
|-----|------|
|scriptName: |Python script file path|
|scriptParamsTemplates: |Script parameters templates array ({0} - number of state vector component)|


---


