using MAT
using Random
using Plots
labels = [1; 2; 3; 12; 18; 17; 25; 15; 24]

plot()
s1 = 300
s2 = 250
ms = 4
using LaTeXStrings
x_label = L"\delta"
y_label = "F1 Scores"
deltas = 10 .^LinRange(0,3,10)
for lab = 1:7
    label = labels[lab]

    if lab < 5
        grownum = 200
        seednum = 10
        color = :blue
    elseif lab < 7
        grownum = 2000
        seednum = 50
        color = :red
    else
        grownum = 10000
        seednum = 200
        color = :green
    end

    outputmat = "Output_VaryDelta/Label_$label"*"_$seednum"*"_$grownum.mat"
    mat = matread(outputmat)
    hl = mat["hl_f1"]
    plot!(deltas, hl, grid = false, #label = LabelNames[label],
     markerstrokewidth = 0, markershape = :circle, linewidth = 2, xaxis = :log10,
     size = (s1,s2), markersize = ms,xguidefont=font(18), legend = false,
     xlabel = x_label, ylabel = y_label, legendfont=font(7))
     #color = color)

end

savefig("Plots/AmazonVaryDelta.pdf")
