"""
Function to perform the complete boosting algorithm
Arguments:
stepno::Int -> Number of Boostingsteps which should be performed
y::String -> Name of the endpoint variable
labels::Array{String,1} -> Names of the potential predictor variables
a::Int -> Number of labels for the covariance matrix to start with
x::Int -> Number of labels which should be additionally called
maxvar::Int -> Maximum number of selected variables
maxvar::Int -> Maximum number of boosting steps
"""
function ds_boosting(stepno::Int, y::String, labels::Array{String,1}, a::Int, x::Int, maxvar::Int=stepno)
	myscratch = Boostscratch(1,
							 Array{Float64,1}(),
							 Array{Float64,1}(),
							 Unibeta(Array{Float64,1}(),
									 Array{String,1}()),
							 Covarmat(Array{Float64,2}(),
							 		  Array{String,1}()),
							 Array{Float64,1}(stepno),
							 Array{String, 1}(stepno),
							 Array{String,1}(),
							 Array{String,1}(),
							 0.1,
							 stepno,
							 1,
							 Array{Float64,1}(),
							 a,
							 x,
							 maxvar,
							 y,
							 labels)

	myscratch.pooledunibeta = calc_unibeta(labels,y)
	myscratch.actualbeta = zeros(length(myscratch.pooledunibeta.unibeta))
	myscratch.actualnom = myscratch.pooledunibeta.unibeta

	

	myscratch.wantedlabels = selectionofcovs(myscratch.pooledunibeta, myscratch.x)
	myscratch.pooledcovarmat = calc_covarmat(myscratch.wantedlabels)
	myscratch.usedlabels = deepcopy(myscratch.wantedlabels)



	while ((myscratch.actualstepno <= myscratch.stepno) &&
			(length(unique(myscratch.selections)) < myscratch.maxvar))
		if(myscratch.actualstepno > 1)
			reboost!(myscratch)
		end
		newselval = boost!(myscratch)
		myscratch.wantedlabels = getselections(myscratch, newselval)
		
		myscratch.actualstepno += 1

		expansion = calc_covarmat(myscratch.wantedlabels)
		merge_covarmats!(myscratch, expansion)

	end

	if(myscratch.actualstepno >= 2)
		println(myscratch.actualstepno," bosstingsteps performed \n recent selected labels: ", myscratch.selections)
	end
	myscratch.actualstepno -= 1
	println("--------------------------------------------------------------------------------------------------------------------------------- \n Algorithm ended successfully.")
	if(myscratch.actualstepno == myscratch.stepno)
		println("Reason for termination: Maximum number of steps reached.")
	else
		println("Reason for termination: Maximum number of variables included.")
	end
	println("Number of boostingsteps total: ", myscratch.actualstepno, "\n Selected variables: ", myscratch.selections, "\n ")
end


"""
Function to perform the complete boosting algorithm
Arguments:
stepno::Int -> Number of Boostingsteps which should be performed
y::String -> Name of the endpoint variable
a::Int -> Number of labels for the covariance matrix to start with
x::Int -> Number of labels which should be additionally called
path_RLibrary::String -> Path to R library where R packages opal, dsBaseClient, dsStatsClient and dsModellingClient are saved
url::Array{String,1} -> Array containing URLs to all DataSHIELD-Server as Strings.
user::String -> String containing the User-Name for login
password::String -> String containin either the password or the name of the private key file.
table::Array{String,1} -> Array containing the table names (either one dimensional if always the same or of the same length as the URL array)
servernames::Array{String,1} -> Array containing the server names (same dimension as url).
check::Bool -> If true: It is checked if the variables are standardized, if false not. The default is true.
ignore::Bool -> Only relevant if check=true. If true: the algorithm will continue even if the variables are not standardizes, if false the algorithm will stop. The default is false.
maxvar::Int -> Maximum number of boosting steps
"""
function ds_boosting(stepno::Int, y::String, a::Int, x::Int, path_RLibrary::String, url::Array{String,1},
	user::String, password::String, table::Array{String,1}, servernames::Array, check::Bool=true, ignore::Bool=false, maxvar::Int=stepno)
	
	ds_start(path_RLibrary, url, user, password,table, servernames, check, ignore)

	labels_pre = R"""
	labels <- ds.colnames('D')[[1]]
	labels <- labels[-which(labels==$y[1])]
	"""
	labels = rcopy(labels_pre)

	ds_boosting(stepno, y, labels, a, x, maxvar)
end
