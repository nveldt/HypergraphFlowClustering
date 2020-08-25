using MAT
using Random
using Statistics
include("../src/HyperLocal.jl")

Add path to wherever matrix is store. Loading take a long time.
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

labels = [1; 2; 3; 12; 18; 17; 25; 15; 24]  # Labels for small nodes
lnum = length(labels)

# See outer parameters
seednum = 10
ntimes = 1
epsilon = 1.0
grownum = 200       # How much to grow seed set by using BestNeighbors
deltas = 10 .^LinRange(0,3,10)
enum = length(deltas)

for lab = 1:9
    label = labels[lab]
    T = findall(x->x ==label,NodeLabels)
    nT = length(T)

    # Different seed set sizes depending on dataset
    if lab < 5
        grownum = 200
        seednum = 10
    elseif lab < 7
        grownum = 2000
        seednum = 50
    else
        grownum = 10000
        seednum = 200
    end

    # For each epsilon we store a different .mat file of outputs
    outputmat = "Output_VaryDelta/Label_$label"*"_$seednum"*"_$grownum.mat"
    println(outputmat)

    # Output from HyperLocal
    hl_pr = zeros(enum)
    hl_re = zeros(enum)
    hl_f1 = zeros(enum)
    hl_time = zeros(enum)
    hl_size = zeros(enum)
    newS = zeros(enum)
    hl_cond = zeros(enum)

    # Generate a new seed set
    p = randperm(nT)
    Rstart = T[p[1:seednum]]
    OneHop = get_immediate_neighbors(H,Ht,Rstart)
    Rmore = BestNeighbors(H,d,Rstart,OneHop,grownum)
    R = union(Rmore,Rstart)

    # Force seed nodes to be contained in output set
    Rs = findall(x->in(x,Rstart),R)

    for e = 1:length(deltas)             # Try each seed set with each delta

        delta = deltas[e]

        # Run HyperLocal
        s = time()
        S, lcond = HyperLocal(H,Ht,order,d,R,epsilon,delta,Rs,true)
        hl_time[e] = time()-s
        condS, volS, cutS = tl_cond(H,S,d,1.0,volA,order)
        pr, re, f1 = PRF(T,S)
        hl_pr[e] = pr
        hl_re[e] = re
        hl_f1[e] = f1
        hl_size[e] = length(S)
        hl_cond[e] = condS
        nS =  length(setdiff(S,R))
        newS[e] = nS

        println("$label ($nT): $f1 \t $nS \t $epsilon")
    end

    matwrite(outputmat, Dict("hl_size"=>hl_size, "newS"=>newS, "hl_time"=>hl_time,
    "hl_pr"=>hl_pr, "hl_re"=>hl_re, "hl_f1"=>hl_f1, "hl_cond"=>hl_cond))

end
