import Composites
using Base.Test

# Step 1: Indentify Material Properties
e1 = [181e9,203e9,38.6e9,76e9]
e2 = [10.3e9,11.2e9,8.27e9,5.5e9]
g12 = [7.17e9,8.4e9,4.14e9,2.3e9]
nu12 = [0.28,0.32,0.26,0.34]

# Step 2: Input Ply Information
theta = [0.0,45.0,-45.0,90.0,-45.0,45.0,0.0]
nply = [5,5,5,10,5,5,5]
tply = [0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001]
matlam = [4,4,4,4,4,4,4]

# Step 3: Calculate z-Coordinates for each Ply Group
z = Composites.getz(tply,nply)
ztest = [-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0]*1e-3
for i = 1:length(z)
  @test isapprox(z[i],ztest[i])
end

# Step 4: Determine the In-plane Stiffness Matrix (Q) for each
#         Unique Material (in principal material coordinates)
q = Composites.getQ(e1,e2,g12,nu12)
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

# Step 5: Transform the In-Plane Stiffness Matrix (Q) for each Ply Group into
#         the Structural (Laminate Coordinate System)
qbar = zeros(Float64,3,3,length(matlam))
for i = 1:length(theta)
  qbar[:,:,i] = Composites.rotQ(q[:,:,matlam[i]],theta[i])
end
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

# Step 6: Determine Laminate Stiffness Matrices (A,B, and D)
A,B,D = Composites.getABD(nply,tply,matlam,theta,e1,e2,g12,nu12)
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

# Step 7: Apply Force and Moment Resultants to Analyze Laminate
force = [25000.0,0.0,0.0]
moment = [1000.0,0.0,0.0]
resultantstress = vcat(force,moment)

# Steps 8-9: Obtain abd Matrix, extract a, b, and d
a,b,d = Composites.compliance(A,B,D)
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

# Step 10: Solve for Laminate Mid-plane strains and curvatures
midplanestrain = a*force+b*moment
curvature = b*force+d*moment
epstest = [2.154e-4,-6.996e-5,5.949e-21]
kappatest = [3.887,-2.332,-0.556]
for i = 1:3
  @test isapprox(signif(midplanestrain[i],4),epstest[i],atol = eps(Float64))
  @test round(curvature[i],3) == kappatest[i]
end

# Step 11: Determine Local Strains in Structural Laminate Coordinates
nlam = length(matlam)
localstrain = zeros(Float64,3,nlam+1)
for k = 1:nlam+1
  localstrain[:,k] = midplanestrain+z[k]*curvature
end

# Step 12: Solve for local stresses in structural laminate coordinates
lowerstress = zeros(Float64,3,nlam)
upperstress = zeros(Float64,3,nlam)
for k = 1:nlam
  lowerstress[:,k] = qbar[:,:,k]*localstrain[:,k]
  upperstress[:,k] = qbar[:,:,k]*localstrain[:,k+1]
end

# Step 13: Transform Local Stresses into Principal Material Coordinates
lowerplystress1 = zeros(Float64,3,nlam)
upperplystress1 = zeros(Float64,3,nlam)
for k = 1:nlam
  lowerplystress1[:,k] = Composites.rotstress(lowerstress[:,k],theta[k])
  upperplystress1[:,k] = Composites.rotstress(upperstress[:,k],theta[k])
end

# Step 11-13: Solve for local stresses
lowerplystrain,upperplystrain = Composites.getplystrain(tply,nply,theta,resultantstress,A,B,D)
lowerplystress = Composites.getplystress(lowerplystrain,q,matlam)
upperplystress = Composites.getplystress(upperplystrain,q,matlam)
