
e1 = fill(148e9,2)
e2 = fill(9.65e9,2)
g12 = fill(4.55e9,2)
nu12 = fill(0.3,2)
theta = [0.0,45.0]
# Single value getQ
qtest = [148.87e9  2.91e9  0.0;
           2.91e9  9.71e9  0.0;
           0.0     0.0     4.55e9]
q = LaminateLayup.getQ(e1[1],e2[1],g12[1],nu12[1],theta[1])
for i = 1:size(qtest,1)
  for j = 1:size(qtest,2)
    @test round(q[i,j]/1e9,2) == qtest[i,j]/1e9
  end
end

# Single value getQ with angle
qbartest = [45.65e9 36.55e9 34.79e9;
            36.55e9 45.65e9 34.79e9;
            34.79e9 34.79e9 38.19e9]
qbar = LaminateLayup.getQ(e1[2],e2[2],g12[2],nu12[2],theta[2])
for i = 1:size(qbartest,1)
  for j = 1:size(qbartest,2)
    @test round(qbar[i,j]/1e9,2) == qbartest[i,j]/1e9
  end
end

# Array getQ
qbarvec = LaminateLayup.getQ(e1,e2,g12,nu12,theta)
for i = 1:size(qtest,1)
  for j = 1:size(qtest,2)
    @test round(qbarvec[i,j,1]/1e9,2) == qtest[i,j]/1e9
  end
end
for i = 1:size(qbartest,1)
  for j = 1:size(qbartest,2)
    @test round(qbarvec[i,j,2]/1e9,2) == qbartest[i,j]/1e9
  end
end

# getABD test
nply = [10,10]
tply = [0.0001,0.0001]
matlam = [1,1]
thetalam = [0.0,45.0]
Atest = [194.52 39.46 34.79;
          39.46 55.36 34.79;
          34.79 34.79 42.74]
Btest = [-51.61 16.82 17.40;
          16.82 17.97 17.40;
          17.40 17.40 16.82]
Dtest = [64.84 13.15 11.60;
         13.15 18.45 11.60;
         11.60 11.60 14.25]

A,B,D = LaminateLayup.getABD(nply,tply,matlam,thetalam,e1,e2,g12,nu12)
for i = 1:size(qbartest,1)
  for j = 1:size(qbartest,2)
    @test round(A[i,j]/1e6,2) == Atest[i,j]
    @test round(B[i,j]/1e3,2) == Btest[i,j]
    @test round(D[i,j],2) == Dtest[i,j]
  end
end

alfatest = [13.44  -4.85  -7.14;
            -4.85  41.81 -21.23;
            -7.14 -21.23  64.95]
betatest = [17.07 -6.01 -11.06;
            -6.01 -5.04 -11.06;
           -11.06 -11.06 -24.05]
deltatest = [40.32 -14.56 -21.41;
            -14.56 125.42 -63.68;
            -21.41 -63.68 194.86]

alfa,beta,delta = LaminateLayup.compliance(A,B,D)
for i = 1:size(qbartest,1)
  for j = 1:size(qbartest,2)
    @test round(alfa[i,j]*1e9,2) == alfatest[i,j]
    @test round(beta[i,j]*1e6,2) == betatest[i,j]
    @test round(delta[i,j]*1e3,2) == deltatest[i,j]
  end
end
