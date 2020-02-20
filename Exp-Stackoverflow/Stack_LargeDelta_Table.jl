using MAT
using Random
include("../src/HyperLocal.jl")
include("../include/FlowSeed.jl")
#
## Read in the matrix, if it isn't read in already
using SparseArrays
s = time()
@time M = matread("../data/processed/stackoverflow_answer_H.mat")
LabelMatrix = M["LabelMatrix"]
LabelNames = M["LabelNames"]
MainLabels = M["MainLabels"]
H = M["H"]
order = round.(Int64,vec(sum(H,dims=2)))
d = vec(sum(H,dims=1))
volA = sum(d)
m,n = size(H)
Ht = sparse(H')
toc = time()-s
println("Done loading things into memory in $toc seconds.")

##
Tconds = Vector{Float64}()
Tsizes = Vector{Int64}()
labels = Vector{Int64}()
for lab = 1:length(MainLabels)
    label = MainLabels[lab]
    T = findnz(LabelMatrix[:,label])[1]
    nT = length(T)
    condT, volT, cutT = tl_cond(H,T,d,1.0,volA,order)
    if condT < .2
        # println("$nT \t $condT \t"*LabelNames[label])
        push!(Tconds,condT)
        push!(labels,label)
        push!(Tsizes,nT)
    end
end

using Statistics
ntimes = 1
epsilon = 1.0
delta = 5000.0
outputmat = "Output_Stack/Set45_$(epsilon)_$(delta).mat"

mat = matread(outputmat)
hl_time = round.(mean(mat["hl_time"],dims = 2),digits = 2)
sl_time = round.(mean(mat["sl_time"],dims = 2),digits = 2)
zl_time = round.(mean(mat["zl_time"],dims = 2),digits = 2)
hl_size= round.(mean(mat["hl_size"],dims = 2),digits = 2)
r = round.(mean(mat["r_f1"],dims = 2),digits = 2)
hl = round.(mean(mat["hl_f1"],dims = 2),digits = 2)
b1 = round.(mean(mat["b1_f1"],dims = 2),digits = 2)
b2 = round.(mean(mat["b2_f1"],dims = 2),digits = 2)
nS = round.(mean(mat["newS"],dims = 2),digits = 2)
sl = round.(mean(mat["sl_f1"],dims = 2),digits = 2)
zl = round.(mean(mat["zl_f1"],dims = 2),digits = 2)

for i = 1:length(labels)
    lab = labels[i]
    println(" $(hl[i]) \t $(sl[i]) \t $(zl[i]) \t $(b1[i]) \t $(b2[i])\t $(r[i]) \t $(hl_time[i]) \t $(sl_time[i]) \t $(zl_time[i])\t $(Tconds[i])"*"\t"*LabelNames[lab])
end
