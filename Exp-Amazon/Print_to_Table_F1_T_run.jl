using MAT
using Statistics

# Labels for Amazon datasets
labels = [1; 2; 3; 12; 18; 17; 25; 15; 24]
names = ["Amazon Fashion", "All Beauty", "Appliances",
"Gift Cards", "Magazine Subscriptions", "Luxury Beauty", "Software",
"Industrial and Scientific","Prime Pantry" ]
sizes = [31; 85; 48; 148; 157; 1581; 802; 5334; 4970]
lnum = length(labels)

# See outer parameters
epsilon = 1.0
s1 = 10
s2 = 50
s3 = 200
g1 = 200
g2 = 2000
g3 = 10000
outputmat = "Output/Smallest_$s1"*"_$g1"*"_$epsilon.mat"

println("\nSeeds = ($s1, $s2, $s3), Grow by ($g1, $g2, $g3), Epsilon = $epsilon")
mat = matread(outputmat)
hl_time = round.(mean(mat["hl_time"],dims = 2),digits = 1)
hl_size= round.(mean(mat["hl_size"],dims = 2),digits = 1)
r = round.(mean(mat["r_f1"],dims = 2),digits = 2)
hl = round.(mean(mat["hl_f1"],dims = 2),digits = 2)
b1 = round.(mean(mat["b1_f1"],dims = 2),digits = 2)
b2 = round.(mean(mat["b2_f1"],dims = 2),digits = 2)
nS = round.(mean(mat["newS"],dims = 2),digits = 2)

for i = 1:size(hl_time,1)
    println(names[i]*" & $(sizes[i]) & $(hl_time[i])& $(hl[i]) & $(b1[i]) & $(b2[i])  & $(r[i]) \\\\")
end

## Get medium clusters
outputmat = "Output/Medium_$s2"*"_$g2"*"_$epsilon.mat"
mat = matread(outputmat)
hl_time = round.(mean(mat["hl_time"],dims = 2),digits = 1)
hl_size= round.(mean(mat["hl_size"],dims = 2),digits = 1)
nS = round.(mean(mat["newS"],dims = 2),digits = 2)
r = round.(mean(mat["r_f1"],dims = 2),digits = 2)
hl = round.(mean(mat["hl_f1"],dims = 2),digits = 2)
b1 = round.(mean(mat["b1_f1"],dims = 2),digits = 2)
b2 = round.(mean(mat["b2_f1"],dims = 2),digits = 2)

for i = 1:size(hl_time,1)
    println(names[i+5]*" & $(sizes[i+5]) & $(hl_time[i])& $(hl[i]) & $(b1[i]) & $(b2[i])  & $(r[i]) \\\\")
end


## Get large clusters
outputmat = "Output/Large_$s3"*"_$g3"*"_$epsilon.mat"
mat = matread(outputmat)
hl_time = round.(mean(mat["hl_time"],dims = 2),digits = 1)
hl_size= round.(Int64,mean(mat["hl_size"],dims = 1))
r = round.(mean(mat["r_f1"],dims = 2),digits = 2)
hl = round.(mean(mat["hl_f1"],dims = 2),digits = 2)
b1 = round.(mean(mat["b1_f1"],dims = 2),digits = 2)
b2 = round.(mean(mat["b2_f1"],dims = 2),digits = 2)
nS = round.(mean(mat["newS"],dims = 2),digits = 2)

for i = 1:size(hl_time,1)
    println(names[i+7]*" & $(sizes[i+7]) & $(hl_time[i]) & $(hl[i]) & $(b1[i]) & $(b2[i])  & $(r[i]) \\\\")
end
