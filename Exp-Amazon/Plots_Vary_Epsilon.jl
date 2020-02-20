using MAT

## Colors
C = [1 0 0;
0 .75 0;
0 0 1;
0 .5 .5;
.5 0 .5;
1 .5 0;
0 0 0;
.75 .75 0;
0.1 0.5 0.6;
0 .75 .75;
.4 .5 .3]

using Plots
# Must be kept in this order, unless you are careful what you do
labels = [1; 2; 3; 12; 18; 17; 25; 15; 24]
names = ["Amazon Fashion", "All Beauty", "Appliances",
"Gift Cards", "Magazine Subsciptions", "Luxury Beauty", "Software",
"Industrial and Scientific","Prime Pantry" ]
sizes = [31; 85; 48; 148; 157; 1581; 802; 5334; 4970]
lnum = length(labels)

# See outer parameters
seednum = 10   # 5
grownum = 300  #[200; 300; 500]
epsis = [10.0 1.0 0.1 0.01 0.001]

plot()
for lab = [1;3;4;5]                 # For each cluster
label = labels[lab]
outputmat = "Output_VaryEps/Label_$label"*"_$seednum"*"_$grownum.mat"

mat = matread(outputmat)
hl_time = mean(mat["hl_time"],dims = 2)
hl = round.(mean(mat["hl_f1"],dims = 2),digits = 2)
b1 = round.(mat["b1_f1"],digits = 2)
b2 = round.(mat["b2_f1"],digits = 2)
r = round.(mat["r_f1"],digits = 2)
re = round.(mat["r_re"],digits = 2)
pr = round.(mat["r_pr"],digits = 2)
nS = mat["newS"]

lw = 3
plot!(epsis',hl,linewidth = lw, grid = false, markershape =:circle,
    legend = :bottomleft, color = RGBA(C[lab,1],C[lab,2],C[lab,3],1),
    xaxis=:log10, label = "label $(label)")

end
savefig("Plots/Smallest_VaryEps_seednum_$seednum"*"_grownum_$grownum.pdf")
