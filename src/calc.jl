"""
	calc_covarmat(wantedlabels)

Function, which calculates the pooled covariance matrix of with respect to the relevant labels (for the first matrix calculation).

Returns an Covarmat object (.covmat=covariance matrix, .lables=correponding variable labels).
The function uses the DataSHIELD function ds_cov.

# Arguments
- `wantedlabels::Array{String,1}`: Array which contains the wantedlabels to give information about which covariances need to be calculated

# Returned objects
- `ds_DistributedBoosting.Covarmat`: .covmat=covariance matrix, .lables=correponding variable labels

# Examples
```julia-repl
julia> calc_covarmat(["X1","X2"])
ds_DistributedBoosting.Covarmat([1.0 0.0492964; 0.0492964 1.0], String["X1", "X2"])
```
"""
function calc_covarmat(wantedlabels::Array{String,1},login::DSLogin)
	# Get the saving place for the covariancmatrix - Initialize as 1 as diag is 1
	covmat = ones(length(wantedlabels),length(wantedlabels))

	# Iteration through the variables (-1 as covariance with itself is variance and thus 1)
	for i = 1 : length(wantedlabels) - 1

		# Get name of first variable
		interim1 = wantedlabels[i]
		# Needs to be rput to be useable in RCall syntax
		@rput interim1

		# Iteration through remaining variables
		for j = i + 1 : length(wantedlabels)
			# Get name of second variable
			interim2 = wantedlabels[j]
			# Needs to be rput to be useable in RCall syntax
			@rput interim2

			#Calculate covariance between var 1 and var 2
			cov_res = try
				R"""
					# Get real variablename, combined with tablename
					cov_x <- paste('D',interim1,sep="$")
					cov_y <- paste('D',interim2,sep="$")

					# Get covariances per cohort using DataSHIELD
					covs <- ds.cov(x=cov_x,y=cov_y)

					# Pooling of covariances per cohort using weighted means, robust estimate (number samples -1)
					res <- ncov <- 0
					for(k in 1:length(covs)){
						res <- res + (unlist(covs[[k]][1]) * (unlist(covs[[k]][2])-1))
						ncov <- ncov + unlist(covs[[k]][2])
					}
					res <- res/(ncov-1)
				"""
			catch
				ds_login(login.url, login.user, login.password, 
					login.table, login.servernames)
				R"""
					# Get real variablename, combined with tablename
					cov_x <- paste('D',interim1,sep="$")
					cov_y <- paste('D',interim2,sep="$")

					# Get covariances per cohort using DataSHIELD
					covs <- ds.cov(x=cov_x,y=cov_y)

					# Pooling of covariances per cohort using weighted means, robust estimate (number samples -1)
					res <- ncov <- 0
					for(k in 1:length(covs)){
						res <- res + (unlist(covs[[k]][1]) * (unlist(covs[[k]][2])-1))
						ncov <- ncov + unlist(covs[[k]][2])
					}
					res <- res/(ncov-1)
				"""
			end
			# Save result
			covmat[i,j] =covmat[j,i] = rcopy(cov_res)
		end
	end
	# Put into ds_DistributedBoosting.pooledcovarmat format
	return Covarmat(covmat, wantedlabels)
end

"""
	calc_covarmat!(myscratch)

Function, which calculates the expansions of the pooled covariance matrix
with respect to the relevant labels saved in myscratch.wantedlabels
and expands the saved covariance matrix in myscratch.

Thus, myscratch.pooledcovarmat is expanded and contains the covariances for myscratch.usedlabels and myscratch.wantedlabels.
myscratch.usedlabels is expanded by myscratch.wantedlabels.
The function uses the DataSHIELD function ds_cov.

# Arguments
myscratch::Boostscratch -> Container with data which are relevant to perform the boostinsteps

# Returned objects
- `ds_DistributedBoosting.pooledcovarmat`: .covmat=covariance matrix, .lables=correponding variable labels

# Examples
```julia-repl
julia> mat1 = calc_covarmat(["X1","X2"])
julia> myscratch = Boostscratch(1, Array{Float64,1}(), Array{Float64,1}(), Unibeta(Array{Float64,1}(), Array{String,1}()), mat1, Array{Float64,1}(10), Array{String,1}(), ["X3","X4"], ["X1","X2"], 0.1, 10, 1, Array{Float64,1}(), 2, 2, 10, "Y", ["X1","X2","X3","X4"])
julia> calc_covarmat!(myscratch)
4-element Array{String,1}:
 "X1"
 "X2"
 "X3"
 "X4"
 julia> myscratch.pooledcovarmat
 ds_DistributedBoosting.Covarmat([1.0 0.0492964 0.0668823 0.188088; 0.0492964 1.0 0.523404 0.262267; 0.0668823 0.523404 1.0 0.528047; 0.188088 0.262267 0.528047 1.0], String["X1", "X2", "X3", "X4"])
```
"""
function calc_covarmat!(myscratch::Boostscratch)
	# Get the covariance of the wantedlabels (symmetrical) - Q4
	covmat = calc_covarmat(myscratch.wantedlabels,myscratch.login)
	# Get the saving place for the covariancmatrices for wantedlabels times usedlabels - Initialize as 1 as diag is 1 (potentially not symmetrical)
	# Two for Q2 and Q3, with dimension interchanged
	covmat2 = ones(length(myscratch.wantedlabels), length(myscratch.usedlabels))
	covmat3 = ones(length(myscratch.usedlabels), length(myscratch.wantedlabels))
	# Q1 is already saved in myscratch, containing covariance of usedlabels

	# Get covariances of usedlabels times wantedlabels
	# Iterate through wantedlabels
	for i = 1 : length(myscratch.wantedlabels)
		# Get name of first wantedlabels
		interim1 = myscratch.wantedlabels[i]
		# Needs to be rput to use in RCall
		@rput interim1

		# Iterate through usedlabels
		for j = 1 : length(myscratch.usedlabels)
			# Get name of first usedlabels
			interim2 = myscratch.usedlabels[j]
			# Needs to be rput to use in RCall
			@rput interim2

			# Calculate covariances using DataSHIELD
			cov_res = try
				R"""
					# Get real variablename, combined with tablename
					cov_x <- paste('D',interim1,sep="$")
					cov_y <- paste('D',interim2,sep="$")

					# Get covariances per cohort using DataSHIELD
					covs <- ds.cov(x=cov_x,y=cov_y)

					# Pooling of covariances per cohort using weighted means, robust estimate (number samples -1)
					res <- ncov <- 0
					for(k in 1:length(covs)){
						res <- res + (unlist(covs[[k]][1]) * (unlist(covs[[k]][2])-1))
						ncov <- ncov + unlist(covs[[k]][2])
					}
					res <- res/(ncov-1)
				"""
			catch
				ds_login(myscratch.login.url, myscratch.login.user, myscratch.login.password, 
					myscratch.login.table, myscratch.login.servernames)
				R"""
					# Get real variablename, combined with tablename
					cov_x <- paste('D',interim1,sep="$")
					cov_y <- paste('D',interim2,sep="$")

					# Get covariances per cohort using DataSHIELD
					covs <- ds.cov(x=cov_x,y=cov_y)

					# Pooling of covariances per cohort using weighted means, robust estimate (number samples -1)
					res <- ncov <- 0
					for(k in 1:length(covs)){
						res <- res + (unlist(covs[[k]][1]) * (unlist(covs[[k]][2])-1))
						ncov <- ncov + unlist(covs[[k]][2])
					}
					res <- res/(ncov-1)
				"""
			end
			# Save result - Interchanged index as interchanged dimensions
			covmat2[i,j] = covmat3[j,i] = rcopy(cov_res)
		end
	end
	# Save in myscratch - wantedlabels added to usedlabels to indicate what covariance values are already called
	myscratch.pooledcovarmat.covarmat = vcat(hcat(myscratch.pooledcovarmat.covarmat,covmat3),hcat(covmat2,covmat.covarmat))
	myscratch.pooledcovarmat.labels = unique(vcat(myscratch.usedlabels,myscratch.wantedlabels))
	myscratch.usedlabels = unique(vcat(myscratch.usedlabels,myscratch.wantedlabels))
end


"""
	calc_unibeta(wantedlabels, y)

Function for the calculation of the pooled univariate estimators.

The function uses the DataSHIELD function ds.glm.

#Arguments:
wantedlabels::Array{String,1} -> Array which contains the wantedlabels to give information about which covariances need to be calculated
y::String -> Array which contains the name of the endpoint variable

# Returned objects
- `ds_DistributedBoosting.pooledcovarmat`: .covmat=covariance matrix, .lables=correponding variable labels

# Examples
```julia-repl
julia> calc_unibeta(["X1","X2"],"Y")
...Details on ds.glm() Execution, e.g. Iterations, deviance, final estimates, per covariate...
ds_DistributedBoosting.Unibeta([2.10027, 0.335083], String["X1", "X2"])
```
"""
function calc_unibeta(wantedlabels::Array{String,1},y::String,login::DSLogin)
	# Intitiate saving place for univariable beta estimates
	unibeta = zeros(length(wantedlabels))

	# Iterate through all variables
	for i = 1 : length(wantedlabels)
		# Get variable name
		label_i = wantedlabels[i]
		# Rput to use in RCall
		@rput label_i
		# Rput name of endpoint to use in RCall
		@rput y

		# Calculate y ~ 1 + label_i estimate
		unib = try
			 R"""
				# Get real variablename, combined with tablename
				vary <- paste('D',y,sep="$")
				varx <- paste('D',label_i,sep="$")

				# Define formula (need to be done seperately to achieve one string)
				myformula <- paste(vary, varx, sep = "~")

				# Get effect estimates - several iterations needed
				interim <- ds.glm(myformula, family = 'gaussian')
				res <- (interim$nsubs / (interim$nsubs - 1)) * interim$coefficients[2,1]
			"""
			catch
			ds_login(login.url, login.user, login.password, 
				login.table, login.servernames)
			R"""
				# Get real variablename, combined with tablename
				vary <- paste('D',y,sep="$")
				varx <- paste('D',label_i,sep="$")

				# Define formula (need to be done seperately to achieve one string)
				myformula <- paste(vary, varx, sep = "~")

				# Get effect estimates - several iterations needed
				interim <- ds.glm(myformula, family = 'gaussian')
				res <- (interim$nsubs / (interim$nsubs - 1)) * interim$coefficients[2,1]
			"""
		end
		# Save result
		unibeta[i] = rcopy(unib)
	end

	# Put into ds_DistributedBoosting.pooledunibeta format
	return Unibeta(unibeta,wantedlabels)
end
