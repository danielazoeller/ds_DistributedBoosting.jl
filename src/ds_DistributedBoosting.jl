module ds_DistributedBoosting

using RCall

export Unibeta, Covarmat, Boostscratch 
export ds_loadPkg, ds_login, ds_check, ds_start, ds_logout
export calc_covarmat, calc_covarmat!, calc_unibeta
export boost!, reboost!
export ds_boosting
export selectionofcovs, getselections, getcovarfromtriangular

"""
	Unibeta(unibeta, labels)

Stucture Unibeta, containing effect estimates and corresponding lables.

#Arguments:
- `unibeta::Array{Float64,1}`: univariable effect estimates
- `labels::Array{String,1}`: labels of univariable effect estimates (should be of same length as unibeta and the same order)
"""
mutable struct Unibeta
	unibeta::Array{Float64,1}
	labels::Array{String,1}
end

"""
	Covarmat(covarmat, labels)

Stucture Covarmat, containing 2-dimensional covariance matrix and corresponding lables.

#Arguments:
- `covarmat::Array{Float64,2}`: 2-dimensional covariance matrix (should be symmetrical)
- `labels::Array{String,1}`: labels of 2-dimensional covairance matrix (should be of same length as one-dim of covarmat and of the same order)
"""
mutable struct Covarmat
	covarmat::Array{Float64,2}
	labels::Array{String,1}
end

"""
	DSLogin(url, user, password, table, servernames)

Stucture DSLogin, containing information to loginto DataShield.

#Arguments:
- `url::Array{String,1}`: Vector with urls to DataSHIELD Servers
- `user::String`: Username
- `password::String`: Password or name of key
- `table::String`: Table name to be loaded
- `servernames::Array{String,1}`: Vector with servernames (personal choice)
"""
mutable struct DSLogin
	url::Array{String,1}
	user::String
	password::String
	table::String
	servernames::Array{String,1}
end

"""
	Boostscratch(actualstepno, actualbeta, actualscore,
					pooledunibeta, pooledcovarmat, actualupdates,
					selections, wantedlabels, usedlabels, nu,
					stepno,	actualsel, actualnom, a,
					x, maxvar, y, labels)

Stucture Boostscratch, containing all relevant information of the ds_DistributedBoosting algorithm.

#Arguments:
- `actualstepno::Int`: Actual stepno in boosting algorithm 
- `actualbeta::Array{Float64,1}`: Actual effect estimate vector (same size as labels)
- `actualscore::Array{Float64,1}`: Actual score function value vector (same size as labels)
- `pooledunibeta::Unibeta`: Pooled univariable effect estiamtes (same size as labels)
- `pooledcovarmat::Covarmat`: Actual pooled covariance matrix (expanded in the course of the algorithm)
- `actualupdates::Array{Float64,1}`: Actual update size (before multiplication with nu)
- `selections::Array{String,1}`: Vector of selected variables (same size as stepno, variables can be in the vector more than once)
- `wantedlabels::Array{String,1}`: Vector containing variable names of the ones for which new covariances need to be called
- `usedlabels::Array{String,1}`: Vector containing variable names of the ones where already covariances are called
- `nu::Float64`: Shrinkage parameter (Experience: 0.1)
- `stepno::Int`: Number of boosting steps to perform maximally
- `actualsel::Int`: Position of last selected variable
- `actualnom::Array{Float64,1}`: Actual denominator for boosting
- `a::Int`: Number of variables for which covariances should be called before the first boosting step (Experience: 20)
- `x::Int`: Number of variables for which covariances should be called additionally to the ones selected by the heuristic (Experience: 20)
- `maxvar::Int`: Number of variables to be selected maximally
- `y::String`: Name of endpoint variable
- `labels::Array{String,1}`: Vector containing names of potential predictors to be considered
"""
mutable struct Boostscratch
	actualstepno::Int
	actualbeta::Array{Float64,1}
	actualscore::Array{Float64,1}
	pooledunibeta::Unibeta
	pooledcovarmat::Covarmat
	actualupdates::Array{Float64,1}
	selections::Array{String,1}
	wantedlabels::Array{String,1}
	usedlabels::Array{String,1}
	nu::Float64
	stepno::Int
	actualsel::Int
	actualnom::Array{Float64,1}
	a::Int
	x::Int
	maxvar::Int
	y::String
	labels::Array{String,1}
	login::DSLogin
end


include("Session.jl")
include("calc.jl")
include("getselections.jl")
include("boost.jl")
include("boosting.jl")

end
