"""
	ds_boosting(stepno, y, labels, a, x, nu=0.1, maxvar=stepno)

Function to perform the complete distributed heuristic boosting algorithm using DataSHIELD based on boost! and reboost!

# Arguments
- `stepno::Int`: Number of Boostingsteps which should be performed
- `y::String`: Name of the endpoint variable
- `labels::Array{String,1}`: Names of the potential predictor variables
- `a::Int`: Number of labels for the covariance matrix to start with
- `x::Int`: Number of labels which should be additionally called
- `nu::Float64`: Shrinkage paramter, needs to be between 0 and 1. Default 0.1.
- `maxvar::Int`: Maximum number of selected variables

# Returned objects
- `ds_DistributedBoosting.Boostscratch`: Boostscratch objects containing the results of the algorithm. 
- `Boostscratch.actualbeta`: Final effect estimate vector (mostly zeros) in order of Boostscratch.labels.
- `Boostscratch.selection`: Names of selected variables.
- `Boostscratch.actualstepno`: Number of performed boosting steps.

# Examples
```julia-repl
julia> myscratch = ds_boosting(10, "Y", ["X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"], 3, 5, 0.1)
...Details on ds.glm() Execution, e.g. Iterations, deviance, final estimates, per covariate...
...Trace of boosting - variables selected per boosting step...
Algorithm ended successfully.
Reason for termination: Maximum number of steps reached.
Number of boostingsteps total: 10
 Histroy of selected variables: String["X3", "X3", "X7", "X10", "X3", "X7", "X10", "X3", "X7", "X10"]
 Selected variables: String["X3", "X7", "X10"] (3 different variables).
ds_DistributedBoosting.Boostscratch(10, [0.0, 0.0, -1.25053, 0.0, 0.0, 0.0, 0.861599, 0.0, 0.0, -0.831607], [4.41114, 0.112281, 6.0799, 0.032712, 0.332016, 0.930665, 5.72913, 3.26847, 1.39211, 6.48666], ds_DistributedBoosting.Unibeta([2.10027, 0.335083, -2.47552, -0.180864, -0.576209, 0.96471, 2.43896, 1.80789, -1.17988, -2.2922], String["X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"]), ds_DistributedBoosting.Covarmat([1.0 0.110723 -0.0383974; 0.110723 1.0 0.178274; -0.0383974 0.178274 1.0], String["X3", "X7", "X10"]), [-0.359872, -0.323885, 0.308955, -0.300352, -0.296071, 0.286693, -0.276565, -0.2707, 0.265951, -0.254689], String["X3", "X3", "X7", "X10", "X3", "X7", "X10", "X3", "X7", "X10"], String[], String["X3", "X7", "X10"], 0.1, 10, 10, [2.10027, 0.335083, -2.47552, -0.180864, -0.576209, 0.96471, 2.43896, 1.80789, -1.17988, -2.2922], 3, 5, 10, "Y", String["X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"])
```
"""
function ds_boosting(stepno::Int, y::String, labels::Array{String,1}, a::Int, x::Int, nu::Float64 = 0.1, maxvar::Int=stepno)
	# Make sure that algorithm performs when number of variables is low
	if(length(labels)<a)
		@warn "Number of potential predictors smaller than a. a is reduced to maximum."
		a = length(labels)
	end
	if(length(labels)<x)
		@warn "Number of potential predictors smaller than x. x is reduced to maximum."
		x = length(labels)
	end

	# Inititate Boostscratch-container
	myscratch = Boostscratch(1,
							 Array{Float64,1}(),
							 Array{Float64,1}(),
							 Unibeta(Array{Float64,1}(),
									 Array{String,1}()),
							 Covarmat(Array{Float64,2}(),
							 		  Array{String,1}()),
							 Array{Float64,1}(stepno),
							 Array{String,1}(),
							 Array{String,1}(),
							 Array{String,1}(),
							 nu,
							 stepno,
							 1,
							 Array{Float64,1}(),
							 a,
							 x,
							 maxvar,
							 y,
							 labels)

	# Calculate pooled univariable effect estimates of all variables in labels
	myscratch.pooledunibeta = calc_unibeta(labels,y)
	# Initiate corresponding vectors
	myscratch.actualbeta = zeros(length(myscratch.pooledunibeta.unibeta))
	myscratch.actualnom = myscratch.pooledunibeta.unibeta

	# Get first wantedlabels for calling covariances 
	myscratch.wantedlabels = selectionofcovs(myscratch.pooledunibeta, myscratch.a)
	@warn (length(myscratch.wantedlabels), " covariances are called, this might take some time.")
	myscratch.pooledcovarmat = calc_covarmat(myscratch.wantedlabels)
	myscratch.usedlabels = deepcopy(myscratch.wantedlabels)

	# Perform first boosting step
	newselval = boost!(myscratch)
	# Get wantedlabels after first boosting step
	myscratch.wantedlabels = getselections(myscratch, newselval)
		
	println(myscratch.actualstepno," bosstingsteps performed \n recent selected labels: ", myscratch.selections)

	myscratch.actualstepno += 1

	stopper = 0
	
	# Perform rest of algorithm: End whne max stepno or max varnumber reached
	while ((myscratch.actualstepno <= myscratch.stepno) && (stopper < myscratch.maxvar))
		# If wantedlables is filled, new covariances need to be called, otherwise the algorithm has all information
		# reboost is necessary to get the actual score vectors after the previous updates for the new variables
		if(length(myscratch.wantedlabels) > 0)
			@warn (length(myscratch.wantedlabels), " new covariances are called, this might take some time.")
			calc_covarmat!(myscratch)
			reboost!(myscratch)
		end
		
		# Perform next boosting step
		newselval = boost!(myscratch)
		# Get new wantedlabesl for covariance calls
		myscratch.wantedlabels = getselections(myscratch, newselval)
				
		stopper = length(unique(myscratch.selections))

		println(myscratch.actualstepno," bosstingsteps performed \n recent selected labels: ", myscratch.selections)
		myscratch.actualstepno += 1
	end

	# Needed as stepno is increased at end of loop
	myscratch.actualstepno -= 1

	# Give back important information
	println("--------------------------------------------------------------------------------------------------------------------------------- \n 
			Algorithm ended successfully.")
	if(myscratch.actualstepno == myscratch.stepno)
		println("Reason for termination: Maximum number of steps reached.")
	else
		println("Reason for termination: Maximum number of variables included.")
	end
	println("Number of boostingsteps total: ", myscratch.actualstepno, "\n Histroy of selected variables: ", myscratch.selections, 
				"\n Selected variables: ", unique(myscratch.selections), " (", length(unique(myscratch.selections)), " different variables).")

	# Return complete results
	return myscratch
end


"""
	ds_boosting(stepno, y, a, x, path_RLibrary, url, user, password, table, servernames, check=true, ignore=false, nu=0.1, maxvar=stepno, labels=Array{String,1}())

Wrapper function to perform the complete distributed heuristic boosting algorithm using DataSHIELD including login and assumption checks.

Calls functions ds_start() and ds_boosting().

# Arguments
- `stepno::Int`: Number of Boostingsteps which should be performed
- `y::String`: Name of the endpoint variable
- `a::Int`: Number of labels for the covariance matrix to start with
- `x::Int`: Number of labels which should be additionally called
- `path_RLibrary::String`: Path to R library where R packages opal, dsBaseClient, dsStatsClient and dsModellingClient are saved
- `url::Array{String,1}`: Array containing URLs to all DataSHIELD-Server as Strings.
- `user::String`:  String containing the User-Name for login
- `password::String`: String containin either the password or the name of the private key file.
- `table::Array{String,1}`: Array containing the table names (either one dimensional if always the same or of the same length as the URL array)
- `servernames::Array{String,1}`: Array containing the server names (same dimension as url).
- `check::Bool=true`: If true: It is checked if the variables are standardized, if false not. The default is true.
- `ignore::Bool=false`: Only relevant if check=true. If true: the algorithm will continue even if the variables are not standardizes, if false the algorithm will stop. The default is false.
- `nu::Float64=0.1`: Shrinkage paramter, needs to be between 0 and 1. Default 0.1.
- `maxvar::Int=stepno`: Maximum number of boosting steps
- `labels::Array{String,1}=Array{String,1}()`: Names of potential predictores to be considered. If empty, all potential variables are used.

# Returned objects
- `ds_DistributedBoosting.Boostscratch`: Containing the information of the results of the distributed boosting algorithm.
- Trace of DataSHIELD (e.g. loaded variables).

# Examples
```julia-repl
julia> result = ds_boosting(30, "Y", 20, 20, "C:/Users/Username/Documents/R/win-library/R-Version",
	   ["https://server_url1:port", "https://server_url2:port"],"user","password",
	   ["Projectname.Tablename"],["Server1","Server2"],false)
```
"""
function ds_boosting(stepno::Int, y::String, a::Int, x::Int, path_RLibrary::String, url::Array{String,1},
	user::String, password::String, table::Array{String,1}, servernames::Array, check::Bool=true, ignore::Bool=false, 
	nu::Float64 =0.1, maxvar::Int=stepno, labels::Array{String,1}=Array{String,1}())
	
	ds_start(path_RLibrary, url, user, password,table, servernames, check, ignore, labels)

	labels_pre = R"""
	labels <- ds.colnames('D')[[1]]
	labels <- labels[-which(labels==$y[1])]
	"""

	if(isempty(labels))
		labels = rcopy(labels_pre)
	end

	myscratch = ds_boosting(stepno, y, labels, a, x, nu, maxvar)

	return myscratch
end
