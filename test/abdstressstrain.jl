# materials
e1 = [181e9,203e9,38.6e9,76e9]
e2 = [10.3e9,11.2e9,8.27e9,5.5e9]
g12 = [7.17e9,8.4e9,4.14e9,2.3e9]
nu12 = [0.28,0.32,0.26,0.34]
rho = zeros(Float64,4)
xt = [1500.0,3500.0,1062.0,1400.0]
xc = [1500.0,1540.0,610.0,235.0]
yt = [40.0,56.0,31.0,12.0]
yc = [246.0,150.0,118.0,53.0]
s = [68.0,98.0,72.0,34.0]
mat = Composites.material(e1,e2,g12,nu12,rho,xt,xc,yt,yc,s)

# laminate
matid = [4,4,4,4,4,4,4]
nply = [5,5,5,10,5,5,5]
tply = [0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001]
theta = [0.0,45.0,-45.0,90.0,-45.0,45.0,0.0]
lam = Composites.laminate(matid,nply,tply,theta)

# getz(tply::Array{Float64,1},nply::Array{Int64,1))
z = Composites.getz(tply,nply)
@testset "getz" begin
  ztest = [-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0]*1e-3
  for i = 1:length(z)
    @test isapprox(z[i],ztest[i])
  end
  # getz(lam::laminate)
  @test Composites.getz(lam) == z
end

# getQ(e1::Array{Float64,1},e2::Array{Float64,1},g12::Array{Float64,1},
#   nu12::Array{Float64,1},theta::Array{Float64,1}=zeros(Float64,length(e1)))
q = Composites.getQ(e1,e2,g12,nu12)
@testset "getQ" begin
  qtest1 = [181.811  2.897 0.0;
             2.897 10.346 0.0;
             0.0    0.0   7.17] #GPa
  for i = 1:3
    for j = 1:3
      @test round(q[i,j,1]/1e9,3) == qtest1[i,j]
    end
  end
  qtest3 = [39.167  2.182 0.0;
             2.182  8.392 0.0;
             0.0    0.0   4.14] #GPa
  for i = 1:3
   for j = 1:3
     @test round(q[i,j,3]/1e9,3) == qtest3[i,j]
   end
  end
  for i = 1:length(mat)
    # getQ(e1::Float64,e2::Float64,g12::Float64,nu12::Float64,theta::Float64=0.0)
    @test Composites.getQ(e1[i],e2[i],g12[i],nu12[i]) == q[:,:,i]
    # getQ(mat::material,theta::Float64)
    @test Composites.getQ(mat[i]) == q[:,:,i]
  end
  # getQ(mat::Array{material,1},theta::Array{Float64,1}
  # getQ(mat::Array{material,1},lam::laminate)
  @test Composites.getQ(mat) == q
end

# rotQ(q::Array{Float64,2},theta::Float64)
qbar = zeros(Float64,3,3,length(matid))
for i = 1:length(theta)
  qbar[:,:,i] = Composites.rotQ(q[:,:,matid[i]],theta[i])
end
@testset "rotQ" begin
  qbartest1 = [76.641  1.886 0.0;
                1.886  5.546 0.0;
                0.0    0.0   2.3] #GPa
  qbartest2 = [23.79  19.19  17.774;
               19.19  23.79  17.774;
               17.774 17.774 19.604] #GPa
  qbartest3 = [23.79   19.19  -17.774;
               19.19   23.79  -17.774;
              -17.774 -17.774  19.604] #GPa
  qbartest4 = [5.546  1.886  0.0;
               1.886  76.641 0.0;
               0.0     0.0   2.3] #GPa
  for i = 1:3
    for j = 1:3
      @test round(qbar[i,j,1]/1e9,3) == qbartest1[i,j]
      @test round(qbar[i,j,2]/1e9,3) == qbartest2[i,j]
      @test round(qbar[i,j,3]/1e9,3) == qbartest3[i,j]
      @test round(qbar[i,j,4]/1e9,3) == qbartest4[i,j]
      @test round(qbar[i,j,5]/1e9,3) == qbartest3[i,j]
      @test round(qbar[i,j,6]/1e9,3) == qbartest2[i,j]
      @test round(qbar[i,j,7]/1e9,3) == qbartest1[i,j]
    end
  end
  # rotQ(q11::Float64,q12::Float64,q22::Float64,q66::Float64,theta::Float64)
  for i = 1:size(qbar,3)
    @test Composites.rotQ(q[1,1,matid[i]],q[1,2,matid[i]],q[2,2,matid[i]],
      q[3,3,matid[i]],theta[i]) == qbar[:,:,i]
  end
end

# getABD(matid::Array{Int64,1},nply::Array{Int64,1},tply::Array{Float64,1},
#   theta::Array{Float64,1},q::Array{Float64,3})
A,B,D = Composites.getABD(matid,nply,tply,theta,q)
@testset "getABD" begin
  Atest = [1.298e8 4.215e7 0.0;
           4.215e7 1.298e8 0.0;
           0.0     0.0     4.381e7] #N
  Btest = zeros(3,3)
  Dtest = [288.317 47.549 17.774;
            47.549 75.033 17.774;
            17.774 17.774 49.759] #N*m
  for i = 1:3
    for j = 1:3
      @test isapprox(signif(A[i,j],4),Atest[i,j])
      @test round(B[i,j],3) == Btest[i,j]
      @test round(D[i,j],3) == Dtest[i,j]
    end
  end
  # getABD(lam::laminate,q::Array{Float64,3})
  @test Composites.getABD(lam,q) == (A,B,D)
end

# compliance(A::Array{Float64,2},B::Array{Float64,2},D::Array{Float64,2})
a,b,d = Composites.compliance(A,B,D)
@testset "compliance" begin
  atest = [8.615e-9 -2.798e-9   0.0;
          -2.798e-9  8.615e-9   0.0;
           0.0       0.0        2.283e-8] #N
  btest = zeros(3,3)
  dtest = [3.887e-3 -2.332e-3 -5.556e-4;
          -2.332e-3  1.596e-2 -4.867e-3;
          -5.556e-4 -4.867e-3  2.203e-2] #N*m
  for i = 1:3
    for j = 1:3
      @test isapprox(signif(a[i,j],4),atest[i,j],atol = eps(Float64))
      @test round(b[i,j],3) == btest[i,j]
      @test isapprox(signif(d[i,j],4),dtest[i,j])
    end
  end
end

# rotstress(stress::Array{Float64,1},theta::Float64)
stress = rand(3)
tsigmatest = [0.5 0.5 1.0;0.5 0.5 -1.0;-0.5 0.5 0.0]*stress
@test rotstress(stress,45) == tsigmatest
# rotstrain(strain::Array{Float64,1},theta::Float64)
strain = rand(3)
tepstest = [0.5 0.5 0.5;0.5 0.5 -0.5;-1.0 1.0 0.0]*strain
@test rotstrain(strain,45) == tepstest

# getplystrain(lam::laminate,resultantstrain::Array{Float64,1})
lowerplystrain,upperplystrain = getplystrain(lam,resultantstrain::Array{Float64,1})
# getplystrain(nply::Array{Int64,1},tply::Array{Float64,1},theta::Array{Float64,1},
#   resultantstrain::Array{Float64,1})
@test (lowerplystrain,upperplystrain) == getplystrain(nply,tply,theta,resultantstrain)
# getplystress(plystrain::Array{Float64,2},q::Array{Float64,3},matid::Array{Int64,1})
lowerplystress = getplystress(lowerplystrain,q,matid)
# getplystress(plystrain::Array{Float64,2},q::Array{Float64,3},lam::laminate) = getplystress(plystrain,q,matid)
getplystress(lowerplystrain,q,lam) == lowerplystress
lowerplystresstest = [-570.672 -54.747 -76.129 80.735 46.629 46.014 456.66;
                        11.223  -9.655  -4.223 -7.519  2.922  7.337 -8.386;
                        2.556   20.799 -13.647 -0.639 7.808 -14.96  -1.917] #MPa
for i = 1:size(lowerplystress,1)
  for j = 1:size(lowerplystress,2)
    round(lowerplystress[i,j]/1e6,3) == lowerplystresstest
  end
end
