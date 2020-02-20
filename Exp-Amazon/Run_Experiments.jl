using MAT
using Random
include("../src/HyperLocal.jl")

## Read in the matrix
s = time()
@time M = matread("../data/AmazonReview5core_H.mat")
NodeLabels = vec(M["NodeLabels"])
NodeNames = M["NodeNames"]
LabelNames = M["LabelNames"]
H = M["H"]
order = round.(Int64,vec(sum(H,dims=2)))

## Extract hyperedges with only one node
e = findall(x->x>1,order)
H = H[e,:]
d = vec(sum(H,dims=1))
order = round.(Int64,vec(sum(H,dims=2)))
volA = sum(d)
m,n = size(H)
Ht = sparse(H')
toc = time()-s
println("Done loading things into memory in $toc seconds.")

## Run several sets of experiments
include("Smallest_Experiments.jl")
include("Medium_Experiments.jl")
include("Largest_Experiments.jl")
