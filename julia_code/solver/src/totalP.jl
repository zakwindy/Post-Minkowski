using DelimitedFiles

file = ARGS[1]
arr = readdlm(file, ' ', Float64, '\n')
G, M, L, T = arr[1,1], arr[1,2], arr[1,3], arr[1,4];
body1 = arr[2,2:end];
body2 = arr[3,2:end];
body3 = arr[4,2:end];

px1, py1, pz1 = body1[4], body1[5], body1[6];
px2, py2, pz2 = body2[4], body2[5], body2[6];
px3, py3, pz3 = body3[4], body3[5], body3[6];

println("Px total is ", px1+px2+px3);
println("Py total is ", py1+py2+py3);
println("Pz total is ", pz1+pz2+pz3);
