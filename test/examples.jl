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
YourStartBuffer = 20
YourRunBuffer = 20
YourCheck = true

YourWarning = true
YourShrink = 0.1
YourMaxvar = 15
YourLabels = ["X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"]

# Perform Boosting for all available variables
result = ds_boosting(YourStepno, YourEndpoint, YourStartBuffer, YourRpath, 
            YourRpath, YourUrls,YourUser,YourPassword,YourTable,
            YourServernames,YourCheck)

# Perform Boosting for selected variables
result_sel = ds_boosting(YourStepno, YourEndpoint, YourStartBuffer, YourRpath, 
            YourRpath, YourUrls,YourUser,YourPassword,YourTable,
            YourServernames,YourCheck,YourWarning, 
            YourShrink, YourMaxvar, YourLabels)
