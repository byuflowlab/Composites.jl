import Composites
using Test

@testset "Textbook ABD" begin
    e1 = fill(148e9, 2)
    e2 = fill(9.65e9, 2)
    g12 = fill(4.55e9, 2)
    nu12 = fill(0.3, 2)
    theta = [0.0, 45.0]

    # Single value getQ
    qtest = [148.87e9  2.91e9  0.0;
               2.91e9  9.71e9  0.0;
               0.0     0.0     4.55e9]
    q = Composites.getQ(e1[1], e2[1], g12[1], nu12[1], theta[1])
    for i = 1:size(qtest, 1)
        for j = 1:size(qtest, 2)
            @test round(q[i,j]/1e9, digits=2) == qtest[i,j]/1e9
        end
    end

    # Single value getQ with angle
    qbartest = [45.65e9 36.55e9 34.79e9;
                36.55e9 45.65e9 34.79e9;
                34.79e9 34.79e9 38.19e9]
    qbar = Composites.getQ(e1[2], e2[2], g12[2], nu12[2], theta[2])
    for i = 1:size(qbartest, 1)
        for j = 1:size(qbartest, 2)
            @test round(qbar[i,j]/1e9, digits=2) == qbartest[i,j]/1e9
        end
    end

    # getABD test
    nply = [10, 10]
    tply = [0.0001, 0.0001]
    matid = [1, 1]
    theta = [0.0, 45.0]
    Atest = [194.52 39.46 34.79;
              39.46 55.36 34.79;
              34.79 34.79 42.74]
    Btest = [-51.61 16.82 17.40;
              16.82 17.97 17.40;
              17.40 17.40 16.82]
    Dtest = [64.84 13.15 11.60;
             13.15 18.45 11.60;
             11.60 11.60 14.25]

    q = Composites.getQ.(e1, e2, g12, nu12)
    A, B, D = Composites.getABD(matid, nply, tply, theta, q)
    for i = 1:size(qbartest, 1)
        for j = 1:size(qbartest, 2)
            @test round(A[i,j]/1e6, digits=2) == Atest[i,j]
            @test round(B[i,j]/1e3, digits=2) == Btest[i,j]
            @test round(D[i,j], digits=2) == Dtest[i,j]
        end
    end

    alfatest = [13.44  -4.85  -7.14;
                -4.85  41.81 -21.23;
                -7.14 -21.23  64.95]
    betatest = [17.07  -6.01 -11.06;
                -6.01  -5.04 -11.06;
               -11.06 -11.06 -24.05]
    deltatest = [40.32 -14.56 -21.41;
                -14.56 125.42 -63.68;
                -21.41 -63.68 194.86]

    alfa, beta, delta = Composites.compliance(A, B, D)
    for i = 1:size(qbartest, 1)
        for j = 1:size(qbartest, 2)
            @test round(alfa[i,j]*1e9, digits=2) == alfatest[i,j]
            @test round(beta[i,j]*1e6, digits=2) == betatest[i,j]
            @test round(delta[i,j]*1e3, digits=2) == deltatest[i,j]
        end
    end
end

### Classical Laminate Theory ###

# materials
e1 = [181e9, 203e9, 38.6e9, 76e9]
e2 = [10.3e9, 11.2e9, 8.27e9, 5.5e9]
g12 = [7.17e9, 8.4e9, 4.14e9, 2.3e9]
nu12 = [0.28, 0.32, 0.26, 0.34]
rho = zeros(Float64, 4)
xt = [1500.0, 3500.0, 1062.0, 1400.0]*1e6
xc = [1500.0, 1540.0, 610.0, 235.0]*1e6
yt = [40.0, 56.0, 31.0, 12.0]*1e6
yc = [246.0, 150.0, 118.0, 53.0]*1e6
s = [68.0, 98.0, 72.0, 34.0]*1e6
t = fill(0.0001, 4)
mat = Composites.Material.(e1, e2, g12, nu12, rho, xt, xc, yt, yc, s, t)

# Laminate
matid = [4, 4, 4, 4, 4, 4, 4]
nply = [5, 5, 5, 10, 5, 5, 5]
tply = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]
theta = [0.0, 45.0, -45.0, 90.0, -45.0, 45.0, 0.0]
lam = Composites.Laminate{Float64}(matid, nply, tply, theta)

# getz(tply::Array{Float64,1},nply::Array{Int64,1))
z = Composites.getz(tply, nply)
@testset "getz" begin
    ztest = [-2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0]*1e-3
    for i = 1:length(z)
        @test isapprox(z[i], ztest[i])
    end
    # getz(lam::Laminate)
    @test Composites.getz(lam) == z
end

# getQ(e1::Array{Float64,1},e2::Array{Float64,1},g12::Array{Float64,1},
#   nu12::Array{Float64,1},theta::Array{Float64,1}=zeros(Float64,length(e1)))
q = Composites.getQ.(e1, e2, g12, nu12)
@testset "getQ" begin
    qtest1 = [181.811 2.897 0.0; 2.897 10.346 0.0; 0.0 0.0 7.17] #GPa
    for i = 1:3
        for j = 1:3
            @test round.(q[1][i,j]/1e9, digits=3) == qtest1[i,j]
        end
    end
    qtest3 = [39.167 2.182 0.0; 2.182  8.392 0.0; 0.0 0.0 4.14] #GPa
    for i = 1:3
        for j = 1:3
            @test round.(q[3][i,j]/1e9, digits=3) == qtest3[i,j]
        end
    end
    for i = 1:length(mat)
        # getQ(e1::Float64,e2::Float64,g12::Float64,nu12::Float64,theta::Float64=0.0)
        @test Composites.getQ(e1[i], e2[i], g12[i], nu12[i]) == q[i][:,:]
        # getQ(mat::material,theta::Float64)
        @test Composites.getQ(mat[i]) == q[i][:,:]
    end
    # getQ(mat::Array{material,1},theta::Array{Float64,1}
    # getQ(mat::Array{material,1},lam::Laminate)
    @test Composites.getQ.(mat) == q
end

# rotQ(q::Array{Float64,2},theta::Float64)
qbar = [Composites.rotQ(q[matid[i]], theta[i]) for i in 1:length(theta)]
@testset "rotQ" begin
    qbartest1 = [76.641 1.886 0.0; 1.886 5.546 0.0; 0.0 0.0 2.3] #GPa
    qbartest2 = [23.79 19.19 17.774; 19.19 23.79 17.774; 17.774 17.774 19.604] #GPa
    qbartest3 = [23.79 19.19 -17.774; 19.19 23.79 -17.774; -17.774 -17.774 19.604] #GPa
    qbartest4 = [5.546 1.886 0.0; 1.886 76.641 0.0; 0.0 0.0 2.3] #GPa
    for i = 1:3
        for j = 1:3
            @test round.(qbar[1][i,j]/1e9, digits=3) == qbartest1[i,j]
            @test round.(qbar[2][i,j]/1e9, digits=3) == qbartest2[i,j]
            @test round.(qbar[3][i,j]/1e9, digits=3) == qbartest3[i,j]
            @test round.(qbar[4][i,j]/1e9, digits=3) == qbartest4[i,j]
            @test round.(qbar[5][i,j]/1e9, digits=3) == qbartest3[i,j]
            @test round.(qbar[6][i,j]/1e9, digits=3) == qbartest2[i,j]
            @test round.(qbar[7][i,j]/1e9, digits=3) == qbartest1[i,j]
        end
    end
    # rotQ(q11::Float64,q12::Float64,q22::Float64,q66::Float64,theta::Float64)
    for i = 1:size(qbar, 3)
        @test Composites.rotQ(q[matid[i]][1,1],q[matid[i]][1,2],q[matid[i]][2,2],
            q[matid[i]][3,3], theta[i]) == qbar[i][:,:]
    end
end

# getABD(matid::Array{Int64,1},nply::Array{Int64,1},tply::Array{Float64,1},
#   theta::Array{Float64,1},q::Array{Float64,3})
A,B,D = Composites.getABD(matid, nply, tply, theta, q)
@testset "getABD" begin
    Atest = [1.298e8 4.215e7 0.0; 4.215e7 1.298e8 0.0; 0.0 0.0 4.381e7] #N
    Btest = zeros(3, 3)
    Dtest = [288.317 47.549 17.774; 47.549 75.033 17.774; 17.774 17.774 49.759] #N*m
    for i = 1:3
        for j = 1:3
            @test isapprox(round(A[i,j], sigdigits=4),Atest[i,j])
            @test round(B[i,j], digits=3) == Btest[i,j]
            @test round(D[i,j], digits=3) == Dtest[i,j]
        end
    end
    # getABD(lam::Laminate,q::Array{Float64,3})
    @test Composites.getABD(lam, q) == (A,B,D)
end

# compliance(A::Array{Float64,2},B::Array{Float64,2},D::Array{Float64,2})
a, b, d = Composites.compliance(A, B, D)
@testset "compliance" begin
    atest = [8.615e-9 -2.798e-9 0.0; -2.798e-9  8.615e-9 0.0;0.0 0.0 2.283e-8] #N
    btest = zeros(3,3)
    dtest = [3.887e-3 -2.332e-3 -5.556e-4; -2.332e-3  1.596e-2 -4.867e-3;
        -5.556e-4 -4.867e-3  2.203e-2] #N*m
    for i = 1:3
        for j = 1:3
            @test isapprox(round(a[i,j],sigdigits=4), atest[i,j], atol = eps(Float64))
            @test round(b[i,j], digits=3) == btest[i,j]
            @test isapprox(round(d[i,j],sigdigits=4), dtest[i,j])
        end
    end
end

@testset "rotstress" begin
    # rotstress(stress::Array{Float64,1},theta::Float64)
    stress = rand(3)
    tsigmatest = [0.5 0.5 1.0;0.5 0.5 -1.0;-0.5 0.5 0.0]*stress
    tsigma = Composites.rotstress(stress, 45.0)
    for i = 1:size(tsigma, 1)
        for j = 1:size(tsigma, 2)
            @test isapprox(tsigma[i,j], tsigmatest[i,j])
        end
    end
end
@testset "rotstrain" begin
    # rotstrain(strain::Array{Float64,1},theta::Float64)
    strain = rand(3)
    tepstest = [0.5 0.5 0.5;0.5 0.5 -0.5;-1.0 1.0 0.0]*strain
    teps = Composites.rotstrain(strain, 45.0)
    for i = 1:size(teps, 1)
        for j = 1:size(teps, 2)
            @test isapprox(teps[i,j], tepstest[i,j])
        end
    end
end

force = [25000.0, 0.0, 0.0]
moment = [1000.0, 0.0, 0.0]
resultantstress = vcat(force, moment)
S = vcat(hcat(A, B), hcat(B', D))
resultantstrain = S\resultantstress

lowerplystrain,upperplystrain = Composites.getplystrain(lam, resultantstrain::Array{Float64,1})
lowerplystress = Composites.getplystress(lowerplystrain, q, matid)
upperplystress = Composites.getplystress(upperplystrain, q, matid)
@testset "getplystrain/getplystress" begin
    # getplystrain(lam::Laminate,resultantstrain::Array{Float64,1})
    # getplystrain(nply::Array{Int64,1},tply::Array{Float64,1},theta::Array{Float64,1},
    #   resultantstrain::Array{Float64,1})
    @test (lowerplystrain, upperplystrain) == Composites.getplystrain(nply, tply,
        theta, resultantstrain)
    # getplystress(plystrain::Array{Float64,2},q::Array{Float64,3},matid::Array{Int64,1})
    # getplystress(plystrain::Array{Float64,2},q::Array{Float64,3},lam::Laminate) = getplystress(plystrain,q,matid)
    @test Composites.getplystress(lowerplystrain, q, lam) == lowerplystress
    lowerplystresstest = [-570.672 -54.747 -76.129 80.735 46.629 46.014 456.66;
                          11.223  -9.655  -4.223 -7.519  2.922  7.337 -8.386;
                          2.556   20.799 -13.647 -0.639 7.808 -14.96  -1.917] #MPa
    for i = 1:3
        for j = 1:length(lowerplystress)
            @test round.(lowerplystress[j][i]/1e6, digits = 3) == lowerplystresstest[i,j]
        end
    end
    # getplystress(plystrain::Array{Float64,2},q::Array{Float64,3},lam::Laminate) = getplystress(plystrain,q,matid)
    @test Composites.getplystress(upperplystrain, q, lam) == upperplystress
    upperplystresstest = [-423.91 -34.595 -35.21  -90.646 87.548  66.167 603.422;
                           8.422  -6.257  -1.841  9.644  5.304  10.736 -11.187;
                           1.917  13.647  -6.496  0.639  14.96 -22.112 -2.556] #MPa
    for i = 1:3
        for j = 1:length(upperplystress)
            @test round(upperplystress[j][i]/1e6, digits=3) == upperplystresstest[i,j]
        end
    end
end

@testset "maxstress" begin
    failmaxstresstest = [
        -0.302793   -0.0247108  -0.0251497  -0.0647469   0.0625344   0.0472619   0.431016;
        1.80387     0.147213    0.149828    0.385726   -0.372545   -0.28156    -2.56775  ;
        0.701835   -0.521387   -0.153435    0.803671    0.441966    0.894661   -0.932254 ;
        -0.158906    0.11805     0.0347399  -0.181963   -0.100068   -0.202565    0.211076 ;
        0.0563767   0.401395   -0.191046    0.0187922   0.439999   -0.650348   -0.0751689;
        -0.0563767  -0.401395    0.191046   -0.0187922  -0.439999    0.650348    0.0751689]
    sfmaxstresstest = [
        -3.30258   -40.4682   -39.7618   -15.4447   15.9912   21.1587     2.3201;
        0.554362    6.79287    6.67431    2.59251  -2.68424  -3.55164   -0.389446;
        1.42484    -1.91796   -6.51743    1.24429   2.26262   1.11774   -1.07267;
        -6.29303     8.47099   28.7853    -5.49561  -9.99323  -4.93669    4.73762;
        17.7378      2.49131   -5.23433   53.2135    2.27273  -1.53764  -13.3034;
        -17.7378     -2.49131    5.23433  -53.2135   -2.27273   1.53764   13.3034]
    failmaxstress,sfmaxstress = Composites.getmatfail(upperplystress, mat,
        lam, "maxstress")
    for i = 1:6
        for j = 1:7
            @test round(failmaxstress[j][i], sigdigits=6) == failmaxstresstest[i,j]
            @test round(sfmaxstress[j][i], sigdigits=6) == sfmaxstresstest[i,j]
        end
    end
end

@testset "tsaiwu" begin
    failtsaiwutest = [2.95173,-0.069493,0.0470993,1.17469,0.260911,1.02618,-1.08208]
    sftsaiwutest = [0.413356,2.93886, 4.85947, 0.873163, 2.02113, 0.983487, 1.90494]
    failtsaiwu, sftsaiwu = Composites.getmatfail(upperplystress, mat, lam, "tsaiwu")
    @test isapprox(failtsaiwu, failtsaiwutest, atol=1e-5)
    @test isapprox(sftsaiwu, sftsaiwutest, atol=1e-5)
end

@testset "hashinrotem" begin
    failhashinrotemtest = [
    -0.407623  -0.0391053  -0.0543778   0.0576675   0.0333063   0.0328674   0.326186;
    2.42839    0.232968    0.323953   -0.343551   -0.198421   -0.195806   -1.94323;
    0.880392  -1.02161    -0.284952   -0.392936    0.112031    0.567472   -0.491518;
    -0.050493   0.407417    0.167466    0.0204784  -0.0557789  -0.212765    0.0282125]
    sfhashinrotemtest = [
    -2.45325   -25.572     -18.3899   17.3408   30.0243   30.4253    3.06574 ;
    0.411795    4.29244     3.08687  -2.91077  -5.0398   -5.10711  -0.514606;
    1.06577     0.989369    1.87333   1.59529   2.98766   1.32748   1.42636 ;
    -4.45025     1.56668     2.44364   6.98798  -4.23414  -2.16795   5.9536  ]
    failhashinrotem, sfhashinrotem = Composites.getmatfail(lowerplystress, mat,
        lam, "hashinrotem")
    for i = 1:4
        for j = 1:7
            @test round(failhashinrotem[j][i], sigdigits=6) == failhashinrotemtest[i,j]
            @test round(sfhashinrotem[j][i], sigdigits=6) == sfhashinrotemtest[i,j]
        end
    end
end
