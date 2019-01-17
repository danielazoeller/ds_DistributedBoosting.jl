"""
	selectionsofcovs(unibeta, numberofcovs, usedlabels=Array{String,1}())

Function to find the numberofcovs labels which seem to be the most relevant and should be called next. 

It finds the x labels with the highest score in the scorestatistic.

# Arguments
- `unibeta::Unibeta`: Contains the scores to find the new labels
- `numberofcovs::Int`: Number of covariances which should be called additionally
- `usedlabels::Array{ASCIIString, 1}=Array{String,1}()`: Array with the labels which already been used (no need to ask for them again). Default is empty when first call.

# Returned objects
- `wantedlabels::Array{String,1}`: Containing the variables names which should be called next

# Examples
```julia-repl
julia> unib = calc_unibeta(["X1","X2","X3","X4","X5"],"Y")
julia> selectionofcovs(unib,2,["X3"])
2-element Array{String,1}:
 "X1"
 "X5"
```
"""
function selectionofcovs(unibeta::Unibeta, numberofcovs::Int, usedlabels::Array{String,1} = Array{String,1}())
	# Get Score
	dummy = deepcopy(unibeta)
	dummy.unibeta = dummy.unibeta.^2

	# If there are already usedlabels, remove them from the potential list of variable names (no need to call them again)
	if(!isempty(usedlabels))
		for i = 1 : length(usedlabels)
			select = findfirst(isequal(usedlabels[i]),dummy.labels)
			deleteat!(dummy.unibeta, select)
			deleteat!(dummy.labels, select)
		end
	end

	# Get the numberofcovs new labels to be called
	wantedlabels = Array{String,1}(undef,numberofcovs)
	for i = 1 : length(wantedlabels)
		if length(dummy.unibeta)==0
			break
		end
		# Select the one with the highest Score
		select = findmax(dummy.unibeta)[2]
		wantedlabels[i] = dummy.labels[select]

		# Delete corresponding variable from list for potential variables from now on
		deleteat!(dummy.unibeta, select)
		deleteat!(dummy.labels, select)
	end
	# Return result
	return wantedlabels
end

""" 
	getselections(myscratch, newselval)
Function to find the labels which seem to be the most relevant and should be called next.

It finds the most relevant labels with respect to the heuristic approach and calls the function selectionofcovs() for the x 
labels which should be called additionally.  

# Arguments
- `myscratch::Boostscratch`: Container with data which are relevant to perform the boostinsteps.
- `newselval::Float64`: value which is used as a criteria for choosing the next labels.

# Returned objects
- `wantedlabels::Array{String,1}`: Containing the variables names which should be called next

# Examples
```julia-repl
julia> lab = ["X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"]
julia> myscratch = Boostscratch(1,Array{Float64,1}(),Array{Float64,1}(),Unibeta(Array{Float64,1}(),Array{String,1}()),Covarmat(Array{Float64}(undef,0,0),Array{String,1}()),Array{Float64,1}(10),Array{String,1}(),Array{String,1}(),Array{String,1}(),0.1,10,1,Array{Float64,1}(),1,1,10,"Y",lab)

julia> myscratch.pooledunibeta = calc_unibeta(lab,"Y")
julia> myscratch.actualbeta = zeros(length(myscratch.pooledunibeta.unibeta))
julia> myscratch.actualnom = myscratch.pooledunibeta.unibeta

julia> myscratch.wantedlabels = selectionofcovs(myscratch.pooledunibeta, myscratch.a)
julia> myscratch.pooledcovarmat = calc_covarmat(myscratch.wantedlabels)
julia> myscratch.usedlabels = deepcopy(myscratch.wantedlabels)

julia> myscratch.wantedlabels = getselections(myscratch, 5.0)
3-element Array{String,1}:
 "X7"
 "X10"
 "X1"
```
"""
function getselections(myscratch::Boostscratch, newselval::Float64)
	# Initiate save place for wantedlabels
	wantedlabels = Array{String,1}()

	# Iterate through actual scorevector
	for i = 1 : length(myscratch.actualnom)
		# If actualscore for variable is bigger than current score and no covariance for the variable is called yet, the variable is selected by the heuristic
		if(myscratch.actualnom[i]^2 >= newselval && 
				typeof(findfirst(isequal(myscratch.pooledunibeta.labels[i]), myscratch.pooledcovarmat.labels)) != Int64)
			wantedlabels = vcat(wantedlabels, myscratch.pooledunibeta.labels[i])
		end
	end

	# If additional variables should be called (myscratch.x) and there are potential variabels (thus some called in step beforehand), call buffer of size x
	if(myscratch.x > 0 && length(wantedlabels) > 0)
		wantedlabels = vcat(wantedlabels, selectionofcovs(Unibeta(myscratch.actualnom, myscratch.pooledunibeta.labels), 
															myscratch.x, vcat(myscratch.pooledcovarmat.labels, wantedlabels)))
	end
	
	# Return wantedlabels
	return wantedlabels
end

""" 
	getcovarfromtriangular(cov, index)
Function to obtain the covariances with the selected variable (one line of covariance matrix).

# Arguments
- `cov::Array{Float64,2}`: Covariance Matrix (saved in Boostscatrch.pooledcovarmat)
- `index::Int`: Selected covariate

# Returned objects
- `wantedlabels::Array{String,1}`: Containing the variables names which should be called next

# Examples
```julia-repl
julia> mat = calc_covarmat(["X1","X2","X3","X4","X5"])
julia> mat.covarmat
5×5 Array{Float64,2}:
 1.0        0.0492964  0.0668823  0.188088  0.170362
 0.0492964  1.0        0.523404   0.262267  0.149141
 0.0668823  0.523404   1.0        0.528047  0.281563
 0.188088   0.262267   0.528047   1.0       0.500581
 0.170362   0.149141   0.281563   0.500581  1.0
 julia> getcovarfromtriangular(mat.covarmat, 2)
 5-element Array{Float64,1}:
 0.0492964
 1.0
 0.523404
 0.262267
 0.149141
```
"""
function getcovarfromtriangular(cov::Array{Float64,2}, index::Int)
	return cov[index,:]
end

