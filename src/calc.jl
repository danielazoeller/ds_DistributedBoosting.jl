"""
Function which calculates the pooled covariance matrix of with respect to the relevant labels (for the first matrix calculation)
Arguments:
wantedlabels::Array{String,1} -> Array which contains the wantedlabels to give information about which covariances need to be calculated
"""
function calc_covarmat(wantedlabels::Array{String,1})
	covmat = ones(length(wantedlabels),length(wantedlabels))
	for i = 1 : length(wantedlabels) - 1
		interim1 = wantedlabels[i]
		for j = i + 1 : length(wantedlabels)
			interim2 = wantedlabels[j]
			cov_res = R"""
			cov_x <- paste('D',$interim1,sep="\$")
			cov_y <- paste('D',$interim2,sep="\$")
			covs <- ds.cov(x=cov_x,y=cov_y)
			res <- ncov <- 0
			for(k in 1:length(covs)){
  				res <- res + (unlist(covs[[k]][1]) * unlist(covs[[k]][2]))
  				ncov <- ncov + unlist(covs[[k]][2])
			}
			res <- res/(ncov-1)
			"""
			covmat[i,j] =covmat[j,i] = rcopy(cov_res)
		end
	end
	return Covarmat(covmat, wantedlabels)
end
"""
Function which calculates the expansions of the pooled covariance matrix with respect to the relevant labels
Arguments:
myscratch::Boostscratch -> Container with data which are relevant to perform the boostinsteps
"""
function calc_covarmat!(myscratch::Boostscratch)
	covmat = calc_covarmat(myscratch.wantedlabels)
	covmat2 = ones(length(myscratch.wantedlabels), length(myscratch.usedlabels))
	covmat3 = ones(length(myscratch.usedlabels), length(myscratch.wantedlabels))
	for i = 1 : length(myscratch.wantedlabels)
		interim1 = myscratch.wantedlabels[i]
		for j = 1 : length(myscratch.usedlabels)
			interim2 = myscratch.usedlabels[j]
			cov_res = R"""
			cov_x <- paste('D',$interim1,sep="\$")
			cov_y <- paste('D',$interim2,sep="\$")
			covs <- ds.cov(x=cov_x,y=cov_y)
			res <- ncov <- 0
			for(k in 1:length(covs)){
				res <- res + (unlist(covs[[k]][1]) * unlist(covs[[k]][2]))
				ncov <- ncov + unlist(covs[[k]][2])
			}
			res <- res/(ncov-1)
			"""
			covmat2[i,j] = covmat3[j,i] = rcopy(cov_res)
		end
	end
	myscratch.pooledcovarmat.covarmat = vcat(hcat(myscratch.pooledcovarmat.covarmat,covmat3),hcat(covmat2,covmat.covarmat))
	myscratch.pooledcovarmat.labels = unique(vcat(myscratch.usedlabels,myscratch.wantedlabels))
	myscratch.usedlabels = unique(vcat(myscratch.usedlabels,myscratch.wantedlabels))
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
			vary <- paste('D',$y,sep="\$")
			varx <- paste('D',$wantedlabels[$i],sep="\$")
			myformula <- paste(vary, varx, sep = "~")
			res <- ds.glm(myformula, family = 'gaussian')\$coefficients[2,1]
		"""
		unibeta[i] = rcopy(unib)
	end
	return Unibeta(unibeta,wantedlabels)
end
