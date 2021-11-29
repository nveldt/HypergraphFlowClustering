using MAT
using Random
using Statistics
using StatsBase
include("../src/HyperLocal.jl")

@time M = matread("../data/stackoverflow_answer_H.mat")
LabelMatrix = M["LabelMatrix"]
LabelNames = M["LabelNames"]
H = M["H"]
order = round.(Int64,vec(sum(H,dims=2)))
m,n = size(H)
@show size(H)

esmall = findall(x->x<50,order)
@time Az = WeightedCliqueExpansion(H[esmall,:], order)

matwrite("ZCE_Stack.mat", Dict("Az"=>Az))

@time As = SimpleCliqueExp(H[esmall,:])

matwrite("SCE_Stack.mat", Dict("As"=>As))
