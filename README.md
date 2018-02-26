# ds_DistributedBoosting.jl

This Julia package implements algorithms for performing heuristic componentwise likelihood-based boosting 
on distributed data under data protection constraints for data distributed via DataSHIELD.
The algorithm is solemnly based on aggregated data and no individual-level data are shared.

The data need to be provided on DataSHIELD server and a login (username / password) needs to be provided.
The most important assumption is that the data are standardized beforehand (mean=0, sd=1). 
Please only provide already standardized data, otherwise the results will not be valid.


## Installation:

The `ds_DistributedBoosting` package is installed by:

    Pkg.clone("https://github.com/danielazoeller/ds_DistributedBoosting.jl.git")


## Overview of functions

The following tables provide an overview of the functions of the package, 
together with a short description.
You can find more detailed descriptions for each function using the Julia help mode 
(entered by typing `?` at the beginning of the Julia command prompt).

### Functions for starting and ending a session

Function name    | Short description
---------------- | -----------------
`ds_loadPkg`     | Load needed R-packages.
`ds_login`       | Login into DataSHIELD server.
`ds_check`       | Check if the loaded variables are standardized.
`ds_start`       | Wrapper function for `ds_loadPKg`, `ds_check`, `ds_start`.
`ds_logout`      | Logout of DataSHIELD server.

### Functions for calculating pooled univariable effects estimates and pooled covariances using DataSHIELD.

Function name    | Short description
---------------- | -----------------
`calc_covarmat`  | Calculate pooled covariance matrix using DataSHIELD.
`calc_covarmat!` | Expand pooled covariance matrix using DataSHIELD.
`calc_unibeta`   | Calculate pooled univariable effect estimate using DataSHIELD.

### Functions for extracting relevant informations for next boosting step

Function name           | Short description
----------------        | -----------------
`selectionsofcovs`      | Get variable names of the ones with the next highest score function.
`getselections`         | Get variables names that should be called next according to heuristic (with buffer).
`getcovarfromtriangular`| Extract covriance of all variables with selected covariate from pooled covariance matrix.

### Functions to perform boosting steps

Function name    | Short description
---------------- | -----------------
`boost!`         | Perform a single boosting step.
`reboost!`       | Re-do boosting steps for newly called covariances.

### Functions to perform complete distributed heuristic boosting using DataSHIELD

Function name    | Short description
---------------- | -----------------
`ds_boosting`    | Wrap all relevant steps - Available as single call or as wrapper with login/logout functions.

## Examples

Prerequisite for running the [example code here](test_part/examples.jl) is that the `ds_DistributedBoosting` package is installed:

    Pkg.clone("https://github.com/danielazoeller/ds_DistributedBoosting.jl.git")

Additionally, the standardized data need to be provided on DataSHIELD server and a login needs to be provided. 
You can use the data stored [here](TestData_deployOnDataSHIELD) containing 500 simulated predictors for 500 individuals on 2 servers (250 each). 
You need to install the R-Packages `opal`, `dsBaseClient`, `dsStatsClient`, and `dsModellingClient` in your R-library.

The example code needs to be adjusted to your own settings.

## Reference

The algorithm is explained and evaluated in:

Zöller, D., Lenz, S., Binder, H. (2018). *Distributed multivariable modeling for signature development under data protection constraints*. BioarXiv.
