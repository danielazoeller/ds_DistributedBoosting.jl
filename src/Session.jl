"""
    ds_loadPkg(path_RLibrary)

Function to load DataSHIELD-Client-Packages needed for Boosting.

The functions uses the package RCall.
The packages need to be installed seperately into the user R-library 
and the corresponding server-side packages on the OPAL-Server 
where DataSHIELD is running.

Needed DataSHIELD-Client-Packages:
opal, dsBaseClient, dsStatsClient, dsModellingClient

# Arguments
- `path_RLibrary`::String: Path to R library where R packages opal, dsBaseClient, dsStatsClient and dsModellingClient are saved.

# Returned objects
- `RCall.RObject{RCall.StrSxp}`: Overview of loaded R-packages.

# Examples
```julia-repl
julia> ds_loadPkg("C:/Users/Username/Documents/R/win-library/R-Version")
RCall.RObject{RCall.StrSxp}
 [1] "dsModellingClient" "dsStatsClient"     "dsBaseClient"
 [4] "fields"            "maps"              "spam"
 [7] "grid"              "dotCall64"         "opal"
[10] "mime"              "rjson"             "RCurl"
[13] "bitops"            "stats"             "graphics"
[16] "grDevices"         "utils"             "datasets"
[19] "methods"           "base"
```
"""
function ds_loadPkg(path_RLibrary::String)
    # using RCall
    R"""
    .libPaths($path_RLibrary)

    library(opal)
    library(dsBaseClient)
    library(dsStatsClient)
    library(dsModellingClient)
    """
end

"""
    ds_login(url, user, password, table, servernames)

Function to login into DataSHIELD-Server.

The function uses the package RCall.
The DataSHIELD function datashield.login() is used.
The function ds_loadPkg() needs to be called first.

# Arguments
- `url::Array{String,1}`: Array containing URLs to all DataSHIELD-Server as Strings.
- `user::String`: String containing the User-Name for login
- `password::String`: String containin either the password or the name of the private key file.
- `table::Array{String,1}`: Array containing the table names (either one dimensional if always the same or of the same length as the URL array)
- `servernames::Array{String,1}`: Array containing the server names (same dimension as url).

# Returned objects
- `RCall.RObject{RCall.StrSxp}`: Overview of loaded DataSHIELD-server.
- Trace of DataSHIELD (e.g. loaded variables).

# Examples
```julia-repl
julia> ds_login(["https://server_url1:port", "https://server_url2:port"],
"user","password",["Projectname.Tablename"],["Server1","Server2"])
RCall.RObject{RCall.VecSxp}
Server1
url: https://server_url1:port
name: Server1
version: x.x.x.-xxxx
username: user
rid: xxxx-xxxx-xxx

Server2
url: https://server_url2:port
name: Server2
version: x.x.x.-xxxx
username: user
rid: xxxx-xxxx-xxx
```
"""
function ds_login(url::Array{String,1},user::String,password::String,table::Array{String,1},servernames::Array{String,1})
    # Compliance with DataSHIELD rules
    if(length(table)>1 && length(table)!=length(url))
        error("table must be either one-dimensional or of the same dimension as url")
    end
    if(length(servernames)!=length(url))
        error("servernames must be of the same dimension as url")
    end

    # Login into DataSHIELD using datashield.login(). All variables are assigned.
    R"""
    my_logindata <- data.frame(server=$servernames, url=$url, user=$user, password=$password, table=$table)

    opal_login <- datashield.login(logins=my_logindata, assign=TRUE)
    """
end

"""
    ds_check(serveranz, tolerance=0.0000000000000005, ignore=false)

Function to check if the variables are standardized.

The function uses the package RCall.
The DataSHIELD functions ds.colnames(), ds.dim(), ds.mean(), ds.var().
Firstly, the mean will be checked (should be 0) 
and only if the check is valid, the sd will be checked (should be 1).
The function ds_login() needs to be called first.

# Arguments
- `serveranz::Int`: Number of used servers.
- `tolerance::Float64=0.0000000005`: Tolerance for comparison of mean and sd.
- `ignore::Bool=false`: If true: the algorithm will continue even if the variables are not standardizes, if false the algorithm will stop. The default is false.
- `labels::Array{String,1}=Array{String,1}()`: Labels of variables that should be checked.

# Returned objects
- No return objects.
- Trace includes overview of already checked variables.

# Examples
```julia-repl
julia> ds_check(2)
```
"""
function ds_check(serveranz::Int, ignore::Bool=false, tolerance::Float64=0.0000000005, labels::Array{String,1}=Array{String,1}())
    if(isempty(labels))
		# Get all names of the assigned table D
        names_cov = rcopy(R"ds.colnames('D')[[1]]")
    else 
        names_cov = deepcopy(labels)
    end
        

    # Check if mean approx 0, loop is stopped if the first is not
    mean_values = R"""
    tolerance <-  $tolerance
    dimension <- ds.dim('D',type="combine")[[1]][2]
    res <- FALSE
    cat("Check of mean: ")
    for (i in 1:dimension) {
        cat($names_cov[i],", ")
        if(any(abs(unlist(ds.mean(paste("D",$names_cov[i],sep="$"), type = "split"))) > tolerance)) {
            res <- TRUE
            break
        }
    }
    res
    """

    if(rcopy(mean_values))
        # mean unqual to 0, thus no standarization was performed
        if(ignore)
            @warn("Values need to be standardized, but mean value is unqual to 0 (variance is not checked and check for mean is stopped). The results of the algorithm might be unvalid.")
        else
            error("Values need to be standardized, but mean value is unqual to 0 (variance is not checked and check for mean is stopped). The results of the algorithm might be unvalid and the evaluation has been stopped.")
        end
    else
        # mean approx. 0, so sd is checked for approx 1, loop is stopped if the first is not
        sd_values = R"""
        tolerance <-  $tolerance
        dimension <- ds.dim('D',type="combine")[[1]][2]
        res <- FALSE
        cat("\n Check of sd: ")
        for (i in 1:dimension) {
            cat($names_cov[i],", ")
            if(any(abs(unlist(ds.var(paste("D",$names_cov[i],sep="$"), type = "split"))-1) > tolerance)) {
                res <- TRUE
                break
            }
        }
        res
        """

        if(rcopy(sd_values))
            # sd unqual to 1, thus no standardization was performed
            if(ignore)
                @warn("Values need to be standardized, but sd value is unqual to 1 (mean is equal to 0 but check of variance is stopped). The results of the algorithm might be unvalid.")
            else
                error("Values need to be standardized, but sd value is unqual to 1 (mean is equal to 0 but check of variance is stopped). The results of the algorithm might be unvalid and the evaluation has been stopped.")
            end
        end
    end
end



"""
    ds_start(path_RLibrary, url, user, password, table, servernames, check=true, ignore=false, labels::Array{String,1}=Array{String,1}())

Function to wrap the complete login process.

Calls the functions ds_loadPKG(), ds_login, and optionally ds_check().

# Arguments
- `path_RLibrary::String`: Path to R library where R packages opal, dsBaseClient, dsStatsClient and dsModellingClient are saved
- `url::Array{String,1}`: Array containing URLs to all DataSHIELD-Server as Strings.
- `user::String`: String containing the User-Name for login
- `password::String`: String containin either the password or the name of the private key file.
- `table::Array{String,1}`: Array containing the table names (either one dimensional if always the same or of the same length as the URL array)
- `servernames::Array{String,1}`: Array containing the server names (same dimension as url).
- `check::Bool=true`: If true: It is checked if the variables are standardized, if false not. The default is true.
- `ignore::Bool=false`: Only relevant if check=true. If true: the algorithm will continue even if the variables are not standardizes, if false the algorithm will stop. The default is false.
- `labels::Array{String,1}=Array{String,1}()`: Labels of variables that should be checked.

# Returned objects
- No objects are returned.
- Trace of DataSHIELD (e.g. loaded variables).

# Examples
```julia-repl
julia> ds_start("C:/Users/Username/Documents/R/win-library/R-Version",["https://server_url1:port", "https://server_url2:port"],
"user","password",["Projectname.Tablename"],["Server1","Server2"],false)
```
"""
function ds_start(path_RLibrary::String, url::Array{String,1},user::String,password::String,table::Array{String,1},servernames::Array{String,1}, check::Bool=true, ignore::Bool=false, tolerance::Float64=0.0000000005, labels::Array{String,1}=Array{String,1}())
    # Load needed R-package
    ds_loadPkg(path_RLibrary)

    # Login to DataSHIELD
    ds_login(url, user, password, table, servernames)

    if(check)
        # Checking of standardization
        @warn("The standardization of the variables will be checked (mean=0 (firstly), sd=1). All variables are checked subsequently and the algorithm stops when at least one variable is not standardized. The algrotihm might take some time.")
        ds_check(length(servernames),ignore,tolerance,labels)
    end
end

"""
    ds_logout()

Function to log out of DataSHIELD, where the connection is made with ds_login().

# Arguments
No arguments.

# Returned objects
- `RCall.RObject{RCall.StrSxp}`: Showing the logout server.

# Examples
```julia-repl
julia> ds_logout()
RCall.RObject{RCall.VecSxp}
Server1
NULL

Server2
NULL
```
"""
function ds_logout()
    R"""
    datashield.logout(opal_login)
    """
end