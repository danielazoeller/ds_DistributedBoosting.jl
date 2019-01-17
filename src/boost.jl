"""
	boost!(myscratch)

Function to perform a single boostingstep,
which is called by the wrapper-function 'boosting'.

# Arguments:
- `myscratch::Boostscratch -> Boostscratch-Data-Container for boosting`

# Returned objects
- `newselval::Float64`: Current score function value
- Update of myscratch (performed one boosting step)

# Examples
```julia-repl
julia> lab = ["X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"]
julia> myscratch = Boostscratch(1,Array{Float64,1}(),Array{Float64,1}(),Unibeta(Array{Float64,1}(),Array{String,1}()),Covarmat(Array{Float64,2}(),Array{String,1}()),Array{Float64,1}(10),Array{String,1}(),Array{String,1}(),Array{String,1}(),0.1,10,1,Array{Float64,1}(),1,1,10,"Y",lab)

julia> myscratch.pooledunibeta = calc_unibeta(lab,"Y")
julia> myscratch.actualbeta = zeros(length(myscratch.pooledunibeta.unibeta))
julia> myscratch.actualnom = myscratch.pooledunibeta.unibeta

julia> myscratch.wantedlabels = selectionofcovs(myscratch.pooledunibeta, myscratch.a)
julia> myscratch.pooledcovarmat = calc_covarmat(myscratch.wantedlabels)
julia> myscratch.usedlabels = deepcopy(myscratch.wantedlabels)

julia> boost!(myscratch)
10.490144453484817
```
"""
function boost!(myscratch::Boostscratch)
	# Get score of all variables (should be devided by variance, but due to standardization this is 1)
	myscratch.actualscore = myscratch.actualnom.^2 # ./varx

	# Get the variable with the highest score - this is the selected one
	myscratch.actualsel = findmax(myscratch.actualscore)[2]

	# Save that this is the currently selected one in selections - vector
	myscratch.selections = vcat(myscratch.selections, myscratch.pooledunibeta.labels[myscratch.actualsel])

	# Calculate the corresponding update
	myscratch.actualupdates[myscratch.actualstepno] =
		myscratch.nu * myscratch.actualnom[myscratch.actualsel]

	# Calculate the update of the corresponding effect estimate
	myscratch.actualbeta[myscratch.actualsel] +=
		myscratch.actualupdates[myscratch.actualstepno]
	
	# Get covariances from all potential variables with selected one from covariance matrix
	covar = getcovarfromtriangular(myscratch.pooledcovarmat.covarmat,
				findfirst(isequal(myscratch.pooledunibeta.labels[myscratch.actualsel]),
					myscratch.pooledcovarmat.labels
				)
			)
	
	# Update corresponding score vector
	@simd for i = 1 : length(covar)
		@inbounds myscratch.actualnom[findfirst(isequal(myscratch.pooledcovarmat.labels[i]),
												myscratch.pooledunibeta.labels)] -=
			myscratch.actualupdates[myscratch.actualstepno] * covar[i]
	end
	# Calculate the current score of selected variable
	newselval = myscratch.actualnom[myscratch.actualsel]^2
	# Return current score of selected variable
	return newselval
end

"""
	reboost!(myscratch)

Function to repeat the boostingsteps (which were already performed)
with new data resulting from a new data call.

# Arguments
- `myscratch::Boostscratch`: Bootscratch-Data-Container for boosting

# Returned objects
- No objects returned.
- `myscratch.actualnom` is updated.

# Examples
```julia-repl
julia> lab = ["X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"]
julia> myscratch = Boostscratch(1,Array{Float64,1}(),Array{Float64,1}(),Unibeta(Array{Float64,1}(),Array{String,1}()),Covarmat(Array{Float64}(undef,0,0),Array{String,1}()),Array{Float64,1}(10),Array{String,1}(),Array{String,1}(),Array{String,1}(),0.1,10,1,Array{Float64,1}(),2,1,10,"Y",lab)

julia> myscratch.pooledunibeta = calc_unibeta(lab,"Y")
julia> myscratch.actualbeta = zeros(length(myscratch.pooledunibeta.unibeta))
julia> myscratch.actualnom = myscratch.pooledunibeta.unibeta

julia> myscratch.wantedlabels = selectionofcovs(myscratch.pooledunibeta, myscratch.a)
julia> myscratch.pooledcovarmat = calc_covarmat(myscratch.wantedlabels)
julia> myscratch.usedlabels = deepcopy(myscratch.wantedlabels)

julia> boost!(myscratch)
julia> myscratch.wantedlabels = ["X9","X10"]
julia> calc_covarmat!(myscratch)
julia> reboost!(myscratch)
"""
function reboost!(myscratch::Boostscratch)
	# Repeat all previous boosting steps
	for i = 1 : (myscratch.actualstepno-1)
		# Get covariance of all variables with the one selected in boosting step i
		covar = getcovarfromtriangular(myscratch.pooledcovarmat.covarmat,
					findfirst(isequal(myscratch.selections[i]), myscratch.pooledcovarmat.labels))
		
		# Update the score functions one the newly called variables
		for j = 1 : length(myscratch.wantedlabels)
			myscratch.actualnom[findfirst(isequal(myscratch.wantedlabels[j]),myscratch.pooledunibeta.labels)] -=
				myscratch.actualupdates[i] * covar[findfirst(isequal(myscratch.wantedlabels[j]), myscratch.pooledcovarmat.labels)]
		end
	end
end
