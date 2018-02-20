"""
Function to find the labels which seem to be the most relevant and should be called next. 
It finds the x labels with the highest score in the scorestatistic.
Arguments:
unibeta::Unibeta -> Contains the scores to find the new labels
numberofcovs::Int -> Number of covariances which should be called additionally
usedlabels::Array{ASCIIString, 1} -> Array with the labels which already been used (no need to ask for them again)
"""
function selectionofcovs(unibeta::Unibeta, numberofcovs::Int, usedlabels::Array{String,1} = Array{String,1}())
	dummy = deepcopy(unibeta)
	dummy.unibeta = dummy.unibeta.^2
	if(!isempty(usedlabels))
		for i = 1 : length(usedlabels)
			select = findfirst(dummy.labels, usedlabels[i])
			deleteat!(dummy.unibeta, select)
			deleteat!(dummy.labels, select)
		end
	end
	wantedlabels = Array{String,1}(numberofcovs)
	for i = 1 : length(wantedlabels)
		select = indmax(dummy.unibeta)
		deleteat!(dummy.unibeta, select)
		wantedlabels[i] = dummy.labels[select]
		deleteat!(dummy.labels, select)
	end
	return wantedlabels
end
""" 
Function to find the labels which seem to be the most relevant and should be called next.
It finds the most relevant labels with respect to the heuristic approach and calls the function selectionofcovs() for the x 
labels which should be called additionally.  
Arguments:
myscratch::Boostscratch -> Container with data which are relevant to perform the boostinsteps
newselval::Float64 -> value which is used as a criteria for choosing the next labels
"""
function getselections(myscratch::Boostscratch, newselval::Float64)
	wantedlabels = Array{String,1}()
	for i = 1 : length(myscratch.actualnom)
		if(myscratch.actualnom[i]^2 >= newselval && findfirst(myscratch.pooledcovarmat.labels, myscratch.pooledunibeta.labels[i]) == 0)
			wantedlabels = vcat(wantedlabels, myscratch.pooledunibeta.labels[i])
		end
	end
	if(myscratch.x > 0 && length(wantedlabels) > 0)
		wantedlabels = vcat(wantedlabels, selectionofcovs(Unibeta(myscratch.actualnom, myscratch.pooledunibeta.labels), 
															myscratch.x, vcat(myscratch.pooledcovarmat.labels, wantedlabels)))
	end
	return wantedlabels
end

""" 
Function to obtain the covariances with the selected variable
Arguments:
cov::Array{Float64,2} -> Covariance Matrix (saved in Boostscatrch.pooledcovarmat)
index::Int -> Selected covariate
"""
function getcovarfromtriangular(cov::Array{Float64,2}, index::Int)
	return cov[index,:]
end

