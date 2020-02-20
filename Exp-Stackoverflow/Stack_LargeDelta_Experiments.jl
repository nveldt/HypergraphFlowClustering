using MAT
using Random
include("../src/HyperLocal.jl")
include("../include/FlowSeed.jl")

## Read in the matrix, if it isn't read in already

s = time()
@time M = matread("../datastackoverflow_answer_H.mat")
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
        println("$nT \t $condT \t"*LabelNames[label])
        push!(Tconds,condT)
        push!(labels,label)
        push!(Tsizes,nT)
    end
end

##
lnum = length(labels)

# See outer parameters
ntimes = 1
epsilon = 1.0
delta = 5000.0
# seednum = 100
# grownum = 10000

# Output from HyperLocal, large delta
hl_pr = zeros(lnum,ntimes)
hl_re = zeros(lnum,ntimes)
hl_f1 = zeros(lnum,ntimes)
hl_time = zeros(lnum,ntimes)
hl_size = zeros(lnum,ntimes)
newS = zeros(lnum,ntimes)
hl_cond = zeros(lnum,ntimes)

# Output from HyperLocal, delta = 1.0
hl_pr1 = zeros(lnum,ntimes)
hl_re1 = zeros(lnum,ntimes)
hl_f11 = zeros(lnum,ntimes)
hl_time1 = zeros(lnum,ntimes)
hl_size1 = zeros(lnum,ntimes)
newS1 = zeros(lnum,ntimes)
hl_cond1 = zeros(lnum,ntimes)

# Output from first baseline 1
b1_pr = zeros(lnum,ntimes)
b1_re = zeros(lnum,ntimes)
b1_f1 = zeros(lnum,ntimes)
b1_cond = zeros(lnum,ntimes)

# Output from baseline 2
b2_pr = zeros(lnum,ntimes)
b2_re = zeros(lnum,ntimes)
b2_f1 = zeros(lnum,ntimes)
b2_cond = zeros(lnum,ntimes)

# Keep track of R
r_pr =  zeros(lnum,ntimes)
r_re =  zeros(lnum,ntimes)
r_f1 =  zeros(lnum,ntimes)
r_cond =  zeros(lnum,ntimes)

sl_pr = zeros(lnum,ntimes)
sl_re = zeros(lnum,ntimes)
sl_f1 =zeros(lnum,ntimes)
sl_cond = zeros(lnum,ntimes)
sl_size = zeros(lnum,ntimes)
sl_time = zeros(lnum,ntimes)

zl_pr = zeros(lnum,ntimes)
zl_re = zeros(lnum,ntimes)
zl_f1 = zeros(lnum,ntimes)
zl_cond = zeros(lnum,ntimes)
zl_size = zeros(lnum,ntimes)
zl_time = zeros(lnum,ntimes)

# For each epsilon we store a different matrix of outputs
outputmat = "Output_Stack/Set45_$(epsilon)_$(delta).mat"
println(outputmat)

maz = matread("ZCE_Stack.mat")
Az = maz["Az"]
mas = matread("SCE_Stack.mat")
As = mas["As"]

## For a fixed epsilon, seednum, and grownum, run experiments
# on each cluster multiple times

for lab = 1:length(labels)
    label = labels[lab]
    T = findnz(LabelMatrix[:,label])[1]
    nT = length(T)
    @show nT,label

    for index = 1:ntimes

        # Generate a new seed set
        seednum = round(Int64,nT/20)
        grownum = round(Int64,min(nT*2,n))
        p = randperm(nT)
        Rstart = T[p[1:seednum]]
        OneHop = get_immediate_neighbors(H,Ht,Rstart)
        Rmore = GrowRpercent(H,d,Rstart,OneHop,grownum)
        R = union(Rmore,Rstart)
        Rs = findall(x->in(x,Rstart),R)     # Force seed nodes to be in output set
        prr, rer, f1r = PRF(T,R)
        r_pr[lab,index] = prr
        r_re[lab,index] = rer
        r_f1[lab,index] = f1r
        condR, volR, cutR = tl_cond(H,R,d,1.0,volA,order)
        r_cond[lab,index] = condR

        # Run HyperLocal with delta = 1.0
        s = time()
        S, lcond = HyperLocal(H,Ht,order,d,R,epsilon,1.0,Rs,true)
        hl_time1[lab,index] = time()-s
        condS, volS, cutS = tl_cond(H,S,d,1.0,volA,order)
        pr, re, f1_d1 = PRF(T,S)
        hl_pr1[lab,index] = pr
        hl_re1[lab,index] = re
        hl_f11[lab,index] = f1_d1
        hl_size1[lab,index] = length(S)
        hl_cond1[lab,index] = condS
        nS =  length(setdiff(S,R))
        newS1[lab,index] = nS

        # Run HyperLocal with delta = 1000
        s = time()
        S, lcond = HyperLocal(H,Ht,order,d,R,epsilon,delta,Rs,true)
        hl_time[lab,index] = time()-s
        condS, volS, cutS = tl_cond(H,S,d,1.0,volA,order)
        pr, re, f1 = PRF(T,S)
        hl_pr[lab,index] = pr
        hl_re[lab,index] = re
        hl_f1[lab,index] = f1
        hl_size[lab,index] = length(S)
        hl_cond[lab,index] = condS
        nS =  length(setdiff(S,R))
        newS[lab,index] = nS

        # First baseline
        kS = nT-length(Rstart)
        B1 = GrowRpercent(H,d,Rstart,OneHop,kS)
        pr1, re1, f11 = PRF(T,B1)
        b1_pr[lab,index] = pr1
        b1_re[lab,index] = re1
        b1_f1[lab,index] = f11

        # Baseline 2
        B2 = GrowR(H,Rstart,OneHop,kS)
        pr2, re2, f12 = PRF(T,B2)
        b2_pr[lab,index] = pr2
        b2_re[lab,index] = re2
        b2_f1[lab,index] = f12

        # Simple Clique Expansion
        nR = length(R)
        Rs_vec = zeros(length(R))
        Rs_vec[Rs] .= 1
        starter = time()
        SL, lcond = FlowSeed(As,R,epsilon,zeros(nR),Rs_vec)
        sl_time[lab,index] = time()-starter
        pr3, re3, f13 = PRF(T,SL)
        condS, volS, cutS = tl_cond(H,SL,d,1.0,volA,order)
        sl_pr[lab,index] = pr3
        sl_re[lab,index] = re3
        sl_f1[lab,index] = f13
        sl_size[lab,index] = length(SL)
        sl_cond[lab,index] = condS

        # Zhou Clique Expansion + SimpleLocal
        starter = time()
        ZL, lcond = FlowSeed(Az,R,epsilon,zeros(nR),Rs_vec)
        zl_time[lab,index] = time()-starter
        pr4, re4, f14 = PRF(T,ZL)
        condS, volS, cutS = tl_cond(H,ZL,d,1.0,volA,order)
        zl_pr[lab,index] = pr4
        zl_re[lab,index] = re4
        zl_f1[lab,index] = f14
        zl_size[lab,index] = length(ZL)
        zl_cond[lab,index] = condS

        println("$label ($nT):\n $f1r \n $f11 \n $f12 \n $f1 \n $f1_d1\n \t $nS "*LabelNames[label])
    end

    matwrite(outputmat, Dict("hl_size"=>hl_size, "newS"=>newS, "hl_time"=>hl_time,
    "hl_pr"=>hl_pr, "hl_re"=>hl_re, "hl_f1"=>hl_f1, "hl_cond"=>hl_cond,
    "hl_pr1"=>hl_pr1, "hl_re1"=>hl_re1, "hl_f11"=>hl_f11, "hl_cond1"=>hl_cond1,
    "sl_pr"=>sl_pr, "sl_re"=>sl_re, "sl_f1"=>sl_f1, "sl_cond"=>sl_cond, "sl_time"=>sl_time,
    "zl_pr"=>zl_pr, "zl_re"=>zl_re, "zl_f1"=>zl_f1, "zl_cond"=>zl_cond, "zl_time"=>zl_time,
    "sl_size"=>sl_size, "zl_size"=>zl_size,
    "b1_pr"=>b1_pr, "b1_re"=>b1_re, "b1_f1"=>b1_f1,
    "r_pr"=>r_pr, "r_re"=>r_re, "r_f1"=>r_f1,
    "b2_pr"=>b2_pr, "b2_re"=>b2_re, "b2_f1"=>b2_f1))

end
