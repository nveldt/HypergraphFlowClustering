using MAT
using Random
include("../include/FlowSeed.jl")
include("../src/HyperLocal.jl")

## Read in the matrix
s = time()
@time M = matread("../data/stackoverflow_answer_H.mat")
LabelMatrix = M["LabelMatrix"]
LabelNames = M["LabelNames"]
MainLabels = M["MainLabels"]
H = M["H"]
order = round.(Int64,vec(sum(H,dims=2)))
d = vec(sum(H,dims=1))
volA = sum(d)
m,n = size(H)
toc = time()-s
println("Done loading things into memory in $toc seconds.")

# Updated MainLabels be the set of topics between 2000-10000
#   (already stored in current MainLabels), to be only those with conductanc < .2
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

##

using Plots
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

ln = LabelNames[labels]
plot()
p = (sortperm(vec(hl)))
y1 = hl[p]
yb1 = b1[p]
yr = r[p]
yb2 = b2[p]
y2 = max.(yb1,yb2)
y4 = max.(sl[p],zl[p])
x = 1:length(y1)
xlabels = ln[p]

## Wide plot
s1 = 750
s2 = 400
ms = 6
stepy = 1
scatter(x,y1,grid = false, markersize = ms, label = "HyperLocal",
markerstrokewidth = 0, markershape = :circle, linewidth = 0,
legend = :topleft, size = (s1, s2), ymirror = false, ylabel = "F1 Scores",
xticks = (1:stepy:length(ln), "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t".*ln[1:stepy:length(ln)]),xrotation = 40,
xtickfont=font(8),
ytickfont=font(11),
guidefont=font(12),
titlefont=font(10),
legendfont=font(9)
)
scatter!(x,y2, markersize = ms,markerstrokewidth = 0,label = "TN/BN",
 markershape = :circle)
 scatter!(x,y4,markersize = ms,markerstrokewidth = 0, label = "FlowSeed",
  markershape = :circle)

savefig("Plots/StackDots_Wide.pdf")
