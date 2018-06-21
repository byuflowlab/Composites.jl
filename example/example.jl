import Composites

# Step 1: Identify Material Properties
e1 = [181e9,203e9,38.6e9,76e9]
e2 = [10.3e9,11.2e9,8.27e9,5.5e9]
g12 = [7.17e9,8.4e9,4.14e9,2.3e9]
nu12 = [0.28,0.32,0.26,0.34]
rho = zeros(Float64,4)
xt = [1500.0,3500.0,1062.0,1400.0]*1e6
xc = [1500.0,1540.0,610.0,235.0]*1e6
yt = [40.0,56.0,31.0,12.0]*1e6
yc = [246.0,150.0,118.0,53.0]*1e6
s = [68.0,98.0,72.0,34.0]*1e6
t = fill(0.0001,4)
mat = Composites.material.(e1,e2,g12,nu12,rho,xt,xc,yt,yc,s,t)

# Step 2: Input Ply Information
matid = [4,4,4,4,4,4,4]
nply = [5,5,5,10,5,5,5]
tply = [0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001]
theta = [0.0,45.0,-45.0,90.0,-45.0,45.0,0.0]
lam = Composites.laminate(matid,nply,tply,theta)

# Step 3: Get Material Stiffness Matrix (Q)
q = Composites.getQ.(mat)

# Step 4: Determine Laminate Stiffness Matrix (A,B, and D)
A,B,D = Composites.getABD(lam,q)

# Step 5: Get Resultant Forces and Moments
force = [25000.0,0.0,0.0]
moment = [1000.0,0.0,0.0]
resultantstress = vcat(force,moment)

# Step 6: Get Resultant Strain
S = vcat(hcat(A,B),hcat(B',D))
resultantstrain = S\resultantstress

# Step 7: Solve for local ply strains and stresses in material coordinate system
lowerplystrain,upperplystrain = Composites.getplystrain(lam,resultantstrain)
lowerplystress = Composites.getplystress(lowerplystrain,q,lam)
upperplystress = Composites.getplystress(upperplystrain,q,lam)

# Step 8: Apply appropriate material failure criteria
lowermaxstress,lsfmaxstress = Composites.getmatfail(lowerplystress,mat,lam,"maxstress")
uppermaxstress,usfmaxstress = Composites.getmatfail(upperplystress,mat,lam,"maxstress")
lowertsai,lsftsai = Composites.getmatfail(lowerplystress,mat,lam,"tsaiwu")
uppertsai,usftsai = Composites.getmatfail(upperplystress,mat,lam,"tsaiwu")
lowerhashinrotem,lsfhashinrotem = Composites.getmatfail(lowerplystress,mat,lam,"hashinrotem")
upperhashinrotem,usfhashinrotem = Composites.getmatfail(upperplystress,mat,lam,"hashinrotem")
