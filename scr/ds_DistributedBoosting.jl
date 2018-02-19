module ds_DistributedBoosting
using RCall

struct Unibeta
	unibeta::Array{Float64,1}
	labels::Array{String,1}
end

struct Covarmat
	covarmat::Array{Float64,2}
	labels::Array{String,1}
end

struct Boostscratch
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
end


include("SessionStart.jl")
include("calc.jl")
include("boost.jl")
include("boosting.jl")


end
