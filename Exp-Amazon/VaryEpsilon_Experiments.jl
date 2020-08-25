labels = [1; 2; 3; 12; 18]
lnum = length(labels)

# See outer parameters
seednum = 10
ntimes = 5
epsis = [10.0 1.0 0.1 0.01 0.001]
delta = 1.0   # stick with the all-or-nothing cut
grownum = 300

enum = length(epsis)

for lab = 2                     # For each cluster
    label = labels[lab]
    T = findall(x->x ==label,NodeLabels)
    nT = length(T)

    # For each epsilon we store a different .mat file of outputs

    outputmat = "Output_VaryEps/Label_$label"*"_$seednum"*"_$grownum.mat"
    println(outputmat)

    # Output from HyperLocal, multiple values of epsilon
    hl_pr = zeros(enum,ntimes)
    hl_re = zeros(enum,ntimes)
    hl_f1 = zeros(enum,ntimes)
    hl_time = zeros(enum,ntimes)
    hl_size = zeros(enum,ntimes)
    newS = zeros(enum,ntimes)
    hl_cond = zeros(enum,ntimes)

    # Output from first baseline 1
    b1_pr = zeros(ntimes)
    b1_re = zeros(ntimes)
    b1_f1 = zeros(ntimes)
    b1_cond = zeros(ntimes)

    # Keep track of R
    r_pr = zeros(ntimes)
    r_re = zeros(ntimes)
    r_f1 = zeros(ntimes)
    r_cond = zeros(ntimes)

    # Output from baseline 2
    b2_pr = zeros(ntimes)
    b2_re = zeros(ntimes)
    b2_f1 = zeros(ntimes)
    b2_cond = zeros(ntimes)

    for index = 1:ntimes                        # Try ntimes different seed sets

        # Generate a new seed set
        p = randperm(nT)
        Rstart = T[p[1:seednum]]
        OneHop = get_immediate_neighbors(H,Ht,Rstart)
        Rmore = BestNeighbors(H,d,Rstart,OneHop,grownum)
        R = union(Rmore,Rstart)
        # Force seed nodes to be in output set
        Rs = findall(x->in(x,Rstart),R)
        prr, rer, f1r = PRF(T,R)
        r_pr[index] = prr
        r_re[index] = rer
        r_f1[index] = f1r
        condR, volR, cutR = tl_cond(H,R,d,1.0,volA,order)
        r_cond[index] = condR

        # First baseline
        kS = nT-length(Rstart)
        B1 = BestNeighbors(H,d,Rstart,OneHop,kS)
        pr1, re1, f11 = PRF(T,B1)
        b1_pr[index] = pr1
        b1_re[index] = re1
        b1_f1[index] = f11

        # Baseline 2
        B2 = TopNeighbors(H,Rstart,OneHop,kS)
        pr2, re2, f12 = PRF(T,B2)
        b2_pr[index] = pr2
        b2_re[index] = re2
        b2_f1[index] = f12

        for e = 1:length(epsis)             # Try each seed set with each epsilon

            epsilon = epsis[e]

            # Run HyperLocal
            s = time()
            S, lcond = HyperLocal(H,Ht,order,d,R,epsilon,delta,Rs,true)
            hl_time[e,index] = time()-s
            condS, volS, cutS = tl_cond(H,S,d,1.0,volA,order)
            pr, re, f1 = PRF(T,S)
            hl_pr[e,index] = pr
            hl_re[e,index] = re
            hl_f1[e,index] = f1
            hl_size[e,index] = length(S)
            hl_cond[e,index] = condS
            nS =  length(setdiff(S,R))
            newS[e,index] = nS

            println("$label ($nT): $f11 \t $f12 \t $f1 \t $nS \t $epsilon")
        end

        matwrite(outputmat, Dict("hl_size"=>hl_size, "newS"=>newS, "hl_time"=>hl_time,
        "hl_pr"=>hl_pr, "hl_re"=>hl_re, "hl_f1"=>hl_f1, "hl_cond"=>hl_cond,
        "b1_pr"=>b1_pr, "b1_re"=>b1_re, "b1_f1"=>b1_f1,"epsis"=>epsis,
        "r_pr"=>r_pr, "r_re"=>r_re, "r_f1"=>r_f1, "r_cond"=>r_cond,
        "b2_pr"=>b2_pr, "b2_re"=>b2_re, "b2_f1"=>b2_f1))

    end
end
