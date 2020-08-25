labels = [15; 24]
lnum = length(labels)

# See outer parameters
delta = 1.0   # all-or-nothing cut
seednum = 200
ntimes = 5
epsis = [1.0]
grownum = 10000

for e = 1:length(epsis)

    # Output from HyperLocal
    hl_pr = zeros(lnum,ntimes)
    hl_re = zeros(lnum,ntimes)
    hl_f1 = zeros(lnum,ntimes)
    hl_time = zeros(lnum,ntimes)
    hl_size = zeros(lnum,ntimes)
    newS = zeros(lnum,ntimes)
    hl_cond = zeros(lnum,ntimes)

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
    r_pr = zeros(lnum,ntimes)
    r_re = zeros(lnum,ntimes)
    r_f1 = zeros(lnum,ntimes)
    r_cond = zeros(lnum,ntimes)

    # For each epsilon we store a different matrix of outputs
    epsilon = epsis[e]
    outputmat = "Output/Large_$seednum"*"_$grownum"*"_$epsilon.mat"
    println(outputmat)

    # For a fixed epsilon, seednum, and grownum, run experiments
    # on each cluster multiple times
    for lab = 1:length(labels)
        label = labels[lab]
        T = findall(x->x ==label,NodeLabels)
        nT = length(T)

        for index = 1:ntimes

            # Generate a new seed set
            p = randperm(nT)
            Rstart = T[p[1:seednum]]
            OneHop = get_immediate_neighbors(H,Ht,Rstart)
            Rmore = BestNeighbors(H,d,Rstart,OneHop,grownum)
            R = union(Rmore,Rstart)
            Rs = findall(x->in(x,Rstart),R)     # Force seed nodes to be in output set
            prr, rer, f1r = PRF(T,R)
            r_pr[index] = prr
            r_re[index] = rer
            r_f1[index] = f1r
            condR, volR, cutR = tl_cond(H,R,d,1.0,volA,order)
            r_cond[index] = condR

            # Run HyperLocal
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
            B1 = BestNeighbors(H,d,Rstart,OneHop,kS)
            pr1, re1, f11 = PRF(T,B1)
            b1_pr[lab,index] = pr1
            b1_re[lab,index] = re1
            b1_f1[lab,index] = f11
            cond, vol, cut = tl_cond(H,B1,d,1.0,volA,order)
            b1_cond[index] = cond


            # Baseline 2
            B2 = TopNeighbors(H,Rstart,OneHop,kS)
            pr2, re2, f12 = PRF(T,B2)
            b2_pr[lab,index] = pr2
            b2_re[lab,index] = re2
            b2_f1[lab,index] = f12
            cond, vol, cut = tl_cond(H,B2,d,1.0,volA,order)
            b2_cond[index] = cond

            println("$label ($nT): $f11 \t $f12 \t $f1 \t $nS")
        end

        matwrite(outputmat, Dict("hl_size"=>hl_size, "newS"=>newS, "hl_time"=>hl_time,
        "hl_pr"=>hl_pr, "hl_re"=>hl_re, "hl_f1"=>hl_f1, "hl_cond"=>hl_cond,
        "b1_pr"=>b1_pr, "b1_re"=>b1_re, "b1_f1"=>b1_f1,"b1_cond"=>b1_cond,"b2_cond"=>b2_cond,
        "r_pr"=>r_pr, "r_re"=>r_re, "r_f1"=>r_f1, "r_cond"=>r_cond,
        "b2_pr"=>b2_pr, "b2_re"=>b2_re, "b2_f1"=>b2_f1))

    end
end
