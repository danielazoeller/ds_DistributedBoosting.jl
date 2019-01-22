using ds_DistributedBoosting

# Install the R-Packages opal, dsBaseClient, dsStatsClient, and dsModellingClient in your R-library

# Please change to your settings!
YourRpath = "C:/Users/Username/Documents/R/win-library/R-Version"
YourUrls = ["https://server_url1:port", "https://server_url2:port"]
YourUser = "user"
YourPassword = "password"
YourTable = ["Projectname.Tablename"]
YourServernames = ["Server1","Server2"]

# Basic parameters
YourStepno = 30
YourEndpoint = "Y"
YourStartBuffer = 2
YourRunBuffer = 2
YourCheck = true
YourTolerance=0.00000000000005

YourWarning = true
YourShrink = 0.1
YourMaxvar = 15
YourLabels = ["X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"]

# Perform Boosting for all available variables
result = ds_boosting(YourStepno, YourEndpoint, YourStartBuffer, YourRunBuffer,
            YourRpath, YourUrls,YourUser,YourPassword,YourTable,
            YourServernames,YourCheck)
# Selected variables are: 
unique(result.selections)
# For the example data, covariances are called 6 times. 
# Of a=2 variables before starting, for additional 5 after 1 boosting step,
# for additional 4 after 7 boosting steps, 4 after 10, 5 after 14 and 3 after 25.
# 9 unique variables are selected:
# X24, X9, X21, X30, X27, X18, X6, X12, X15

# Perform Boosting for selected variables (YourLabels). Check will be performed only for these.
result_sel = ds_boosting(YourStepno, YourEndpoint, YourStartBuffer, YourRunBuffer,
            YourRpath, YourUrls,YourUser,YourPassword,YourTable,
            YourServernames,YourCheck,YourWarning, YourTolerance,
            YourShrink, YourMaxvar, YourLabels)
# For the example data, covariances are called 3 times.
# Of a=2 variables before starting, for additional 4 after 4 boosting steps,
# and 3 after 6 boosting steps.
# 6 unique variables are selected:
# X9, X6, X3, X5, X10, X8