"""
Function to load DataSHIELD-Client-Packages.
The packages need to be installed seperately into the user R-library.
Needed DataSHIELD-Client-Packages:
opal, dsBaseClient, dsStatsClient, dsModellingClient
Arguments:
path_RLibrary::String -> Path to R library where R packages opal, dsBaseClient, dsStatsClient and dsModellingClient are saved
"""
function ds_loadPkg(path_RLibrary::String)
    R"""
    .libPaths($path_RLibrary)

    library(opal)
    library(dsBaseClient)
    library(dsStatsClient)
    library(dsModellingClient)
    """
end

"""
Function to login into DataSHIELD-Server.
Arguments:
url::Array{String,1} -> Array containing URLs to all DataSHIELD-Server as Strings.
user::String -> String containing the User-Name for login
password::String -> String containin either the password or the name of the private key file.
table::Array{String,1} -> Array containing the table names (either one dimensional if always the same or of the same length as the URL array)
servernames::Array{String,1} -> Array containing the server names (same dimension as url).
"""
function ds_login(url::Array{String,1},user::String,password::String,table::Array{String,1},servernames::Array{String,1})
    if(length(table)>1 && length(table)!=length(url))
        error("table must be either one-dimensional or of the same dimension as url")
    end
    if(length(servernames)!=length(url))
        error("servernames must be of the same dimension as url")
    end

    R"""


    my_logindata <- data.frame(server=$servernames, url=$url, user=$user, password=$password, table=$table)

    opal_login <- datashield.login(logins=my_logindata, assign=TRUE)
    """
end

"""
Function to check if the variables are standardized.
Arguments:
serveranz::Int -> Number of used servers.
ignore::Bool -> If true: the algorithm will continue even if the variables are not standardizes, if false the algorithm will stop. The default is false.
"""
function ds_check(serveranz::Int, ignore::Bool=false)
    names_cov = rcopy(R"ds.colnames('D')[[1]]")
    mean_values = R"""
    tolerance <-  .Machine[which(names(.Machine)=="double.eps")]
    dimension <- ds.dim('D',type="combine")[[1]][2]
    res <- FALSE
    cat("Check of mean: ")
    for (i in 1:dimension) {
        cat(i,", ")
        if(any(unlist(ds.mean(paste("D",$names_cov[i],sep="$"), type = "split")) > tolerance)) {
            res <- TRUE
            break
        }
    }
    res
    """
    if(rcopy(mean_values))
        if(ignore)
            warn("Values need to be standardized, but mean value is unqual to 0 (variance is not be checked). The results of the algorithm might be unvalid.")
        else
            error("Values need to be standardized, but mean value is unqual to 0 (variance is not be checked). The results of the algorithm might be unvalid and the evaluation has been stopped.")
        end
    else
        sd_values = R"""
        tolerance <-  .Machine[which(names(.Machine)=="double.eps")]
        dimension <- ds.dim('D',type="combine")[[1]][2]
        res <- FALSE
        cat("Check of sd: ")
        for (i in 1:dimension) {
            cat(i,", ")
            if(any(unlist(ds.var(paste("D",$names_cov[i],sep="$"), type = "split")) > tolerance)) {
                res <- TRUE
                break
            }
        }
        res
        """

        if(rcopy(sd_values))
            if(ignore)
                warn("Values need to be standardized, but sd value is unqual to 1 (mean is equal to 0). The results of the algorithm might be unvalid.")
            else
                error("Values need to be standardized, but sd value is unqual to 1 (mean is equal to 0). The results of the algorithm might be unvalid and the evaluation has been stopped.")
            end
        end
    end
end



"""
Function to wrap the complete login process.
Arguments:
path_RLibrary::String -> Path to R library where R packages opal, dsBaseClient, dsStatsClient and dsModellingClient are saved
url::Array{String,1} -> Array containing URLs to all DataSHIELD-Server as Strings.
user::String -> String containing the User-Name for login
password::String -> String containin either the password or the name of the private key file.
table::Array{String,1} -> Array containing the table names (either one dimensional if always the same or of the same length as the URL array)
servernames::Array{String,1} -> Array containing the server names (same dimension as url).
check::Bool -> If true: It is checked if the variables are standardized, if false not. The default is true.
ignore::Bool -> Only relevant if check=true. If true: the algorithm will continue even if the variables are not standardizes, if false the algorithm will stop. The default is false.
"""
function ds_start(path_RLibrary::String, url::Array{String,1},user::String,password::String,table::Array{String,1},servernames::Array{String,1}, check::Bool=true, ignore::Bool=false)
    ds_loadPkg(path_RLibrary)

    ds_login(url, user, password, table, servernames)
    if(check)
        warn("The standardization of the variables will be checked (mean=0 (firstly), sd=1). All variables are checked subsequently and the algorithm stops when at least one variable is not standardized. The algrotihm might take some time.")
        ds_check(length(servernames),ignore)
    end
end

"""
Function to log out of DataSHIELD.
"""
function ds_logout()
    R"""
    datashield.logout(opal_login)
    """
end