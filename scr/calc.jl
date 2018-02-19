"""
Function which calculates the pooled covariance matrix of with respect to the relevant labels (for the first matrix calculation)
Arguments:
wantedlabels::Array{String,1} -> Array which contains the wantedlabels to give information about which covariances need to be calculated
"""
function calc_covarmat(wantedlabels::Array{String,1})
	covmat = triu(ones(length(wantedlabels),length(wantedlabels)))
	for i = 1 : length(wantedlabels) - 1
		interim1 = wantedlabels[i]
		for j = i + 1 : length(wantedlabels)
			interim2 = wantedlabels[j]
			cov_res = R"""
			cov_x <- paste('D',$interim1,sep="$")
			cov_y <- paste('D',$interim2,sep="$")
			covs <- ds.cov(x=cov_x,y=cov_y)
			res <- ncov <- 0
			for(i in 1:length(covs)){
  				res <- res + (unlist(covs[[i]][1]) * unlist(covs[[i]][2]))
  				ncov <- ncov + unlist(covs[[i]][2])
			}
			res <- res/(ncov-1)
			"""
			covmat[i,j] =covmat[j,i] = rcopy(Array{Float64},cov_res)
		end
		covmat[i,i] = 1.0
	end
	return Covarmat(covmat, wantedlabels)
end
"""
Function which calculates the expansions of the pooled covariance matrix with respect to the relevant labels
Arguments:
wantedlabels::Array{String,1} -> Array which contains the wantedlabels to give information about which covariances need to be calculated
usedlabels::Array{String,1} -> Array which contains the labels which were already used to give information about which covariances need to be calculated
"""
function calc_covarmat(wantedlabels::Array{String,1}, usedlabels::Array{String,1})
	covmat = triu(ones(length(wantedlabels), length(wantedlabels)))
	covmat2 = zeros(length(usedlabels), length(wantedlabels))
	for i = 1 : max(length(wantedlabels), length(labelmap2))
		interim1 = wantedlabels[i]
		for j = i+1 : max(length(wantedlabels), length(usedlabels))
			interim2 = wantedlabels[j]
			cov_res = R"""
			cov_x <- paste('D',$interim1,sep="$")
			cov_y <- paste('D',$interim2,sep="$")
			covs <- ds.cov(x=cov_x,y=cov_y)
			res <- ncov <- 0
			for(i in 1:length(covs)){
				res <- res + (unlist(covs[[i]][1]) * unlist(covs[[i]][2]))
				ncov <- ncov + unlist(covs[[i]][2])
			}
			res <- res/(ncov-1)
			"""
			if(i <= length(usedlabels) && j <= length(wantedlabels))
				@inbounds covmat2[i,j] = covmat2[j,i] = rcopy(Array{Float64},cov_res)
			end
			if(j >= i+1 && i <= length(wantedlabels) && j <= length(wantedlabels))
				@inbounds covmat[i,j]  = covmat2[j,i] = rcopy(Array{Float64},cov_res)
			end
		end
		covmat[i,i] = covmat2[i,i] = 1.0
	end
	return Covarmat(vcat(covmat2, covmat), vcat(usedlabels, wantedlabels))
end
"""
Function for the calculation of the pooled univariate estimators
Arguments:
wantedlabels::Array{String,1} -> Array which contains the wantedlabels to give information about which covariances need to be calculated
y::String -> Array which contains the name of the endpoint variable
"""
function calc_unibeta(wantedlabels::Array{String,1},y::String)
	unibeta = Array{Float64}(length(wantedlabels))
	for i = 1 : length(wantedlabels)
		unib = R"""
			vary <- paste('D',$y,sep="$")
			varx <- paste('D',$wantedlabels[$i],sep="$")
			myformula <- paste(vary, varx, sep = "~")
			res <- ds.glm(myformula, family = 'gaussian')$coefficients[2,1]
		"""
		unibeta[i] = rcopy(unib)
	end
	return Unibeta(unibeta,wantedlabels)
end
