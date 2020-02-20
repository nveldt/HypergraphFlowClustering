using MAT
deltas = 10 .^LinRange(0,3.7,10)
seednum = 100
epsilon = 1.0
grownum = 10000

# M = matread("../data/processed/stackoverflow_answer_H.mat")
LabelMatrix = M["LabelMatrix"]
LabelNames = M["LabelNames"]
H = M["H"]
order = round.(Int64,vec(sum(H,dims=2)))
d = vec(sum(H,dims=1))
volA = sum(d)
m,n = size(H)


labels = [25849; 27596;28918;29386;43507]
lnum = length(labels)
v = [22943;28886;5713]


##
using Plots
plot()
s1 = 300
s2 = 250
ms = 4
lw = 2
using LaTeXStrings
x_label = L"\delta"
y_label = "F1 Scores"
for i = 1:3
    label = labels[i]
    outputmat = "Output_VaryDelta/Cluster_$(label)_$grownum"*"_$seednum"*"_$epsilon.mat"
    mat = matread(outputmat)
    hl_f1 = mat["hl_f1"]
    deltas = mat["deltas"]

    plot!(deltas, hl_f1[i,:],grid = false,label = LabelNames[label],
     markerstrokewidth = 0, markershape = :circle, linewidth = 2, xaxis = :log10,
     size = (s1,s2), markersize = ms,
     xlabel = x_label, ylabel = y_label, legend = :bottomright)
end

## These were run and stored and slightly differently
v = [22943;28886;5713]
for i = 1:3
    label = v[i]
    outputmat = "Output_VaryDelta/Cluster_$(label)_$grownum"*"_$seednum"*"_$epsilon.mat"
    mat = matread(outputmat)
    hl_f1 = mat["hl_f1"]

    plot!(deltas, vec(hl_f1[:,1]),grid = false,label = LabelNames[label],
     markerstrokewidth = 0, markershape = :circle, linewidth = 2, xaxis = :log10,
     size = (s1,s2), markersize = ms,xguidefont=font(18),legend = false,
     xlabel = x_label, ylabel = y_label, legendfont=font(7))

end

savefig("Plots/Stack_VaryDelta.pdf")
