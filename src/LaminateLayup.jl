"""
This module uses classical laminate theory to calculate stresses and strains in
a thin plate laminate.  Analysis is based on chapter 3 of "Design and Analysis
of Composite Structures" by Kassapoglou.
"""
module LaminateLayup
export stiffness,compliance,plystrain,plystress

"""
    stiffness(e1,e2,g12,anu12,n_lamina,n_plies,t_lam,tht_lam,mat_lam)

Returns stiffness matrix for a layup using classical laminated plate theory.

# Arguments:
Material Properties
* `e1::Array{Float64,1}`: E1
* `e2::Array{Float64,1}`: E2
* `g12::Array{Float64,1}`: G12
* `anu12::Array{Float64,1}`: Nu12
Layup Properties
* `n_lamina`: number of lamina
* `n_plies::Array{Int64,1}`: number of plies for each lamina
* `t_lam::Array{Float64,1}`: thickness of a ply (m) for the lamina
* `tht_lam::Array{Float64,1}`: orientation (deg) for each lamina
* `mat_lamU::Array{Int64,1}`: material id for each lamina
"""
function stiffness(e1,e2,g12,anu12,n_lamina,n_plies,t_lam,tht_lam,mat_lam)

  # Implicitly assigned inputs
  if !(length(e1) == length(e2) == length(g12) == length(anu12))
    error("lengths of specified material properties do not match")
  end
  n_materials = length(e1)

  if !(length(n_plies) == length(mat_lam) == length(t_lam) == length(tht_lam))
    error("lengths of lamina properties do not match")
  end

  # material properties
  anud = 1.0 - anu12.*anu12.*e2./e1
  q11 = e1 ./ anud
  q22 = e2 ./ anud
  q12 = anu12.*e2 ./ anud
  q66 = g12

  # Get z-values
  t_lam = t_lam.*n_plies
  z0 = sum(t_lam)./2.0
  z = zeros(n_lamina+1)
  for k = 1:n_lamina
    z[k+1] = sum(t_lam[1:k])
  end
  z = z-z0

  # Loop through layers
  A = zeros(3,3)
  B = zeros(3,3)
  D = zeros(3,3)
  for k = 1:n_lamina
    # Rotate material stiffness properties by specifed angle
    theta = tht_lam[k]
    qxx = q11[mat_lam[k]]
    qyy = q22[mat_lam[k]]
    qxy = q12[mat_lam[k]]
    qss = q66[mat_lam[k]]
    qt = rotmatprop(theta,qxx,qyy,qxy,qss)

    # Find lumped stiffness properties for the laminate
    for i = 1:3
      for j = 1:3
        A[i,j] = A[i,j]+qt[i,j]*(z[k+1]-z[k])
        B[i,j] = B[i,j]+qt[i,j]/2.0*(z[k+1]^2.0-z[k]^2.0)
        D[i,j] = D[i,j]+qt[i,j]/3.0*(z[k+1]^3.0-z[k]^3.0)
      end
    end
  end

  # Assemble and return stiffness matrix
  S = [A B; B' D]
  return S
end #stiffness

"""
    compliance(e1,e2,g12,anu12,n_lamina,n_plies,t_lam,tht_lam,mat_lam)

Returns compliance matrix for a layup using classical laminated plate theory.

# Arguments:
Material Properties
* `e1::Array{Float64,1}`: E1
* `e2::Array{Float64,1}`: E2
* `g12::Array{Float64,1}`: G12
* `anu12::Array{Float64,1}`: Nu12
Layup Properties
* `n_lamina`: number of lamina
* `n_plies::Array{Int64,1}`: number of plies for each lamina
* `t_lam::Array{Float64,1}`: thickness of a ply (m) for the lamina
* `tht_lam::Array{Float64,1}`: orientation (deg) for each lamina
* `mat_lamU::Array{Int64,1}`: material id for each lamina

"""
function compliance(e1::Array{Float64,1}, e2::Array{Float64,1},
  g12::Array{Float64,1}, anu12::Array{Float64,1},
  n_lamina::Int64,n_plies::Array{Int64,1}, t_lam::Array{Float64,1},
  tht_lam::Array{Float64,1}, mat_lam::Array{Int64,1})

  S = stiffness(e1,e2,g12,anu12,n_lamina,n_plies,t_lam,tht_lam,mat_lam)
  A = S[1:3,1:3]
  invA = inv(A)
  delta = (D-B*invA*B)^-1
  beta = -A*B*delta
  alfa = invA+invA*B*delta*B*invA
  C = [alfa,beta,beta',delta]
  return C
end # compliance

"""
    plystrain(S,Nx,Ny,Nxy,Mx,My,Mxy,n_lamina,n_plies,t_lam,tht_lam)

Calculate strains in each ply aligned with principal material direction
Returns epsxU,epsxL,epsyU,epsyL,gammaxyU,gammaxyL
"""
function plystrain(S,Nx,Ny,Nxy,Mx,My,Mxy,
  n_lamina::Int64,n_plies::Array{Int64,1}, t_lam::Array{Float64,1}, tht_lam)

  # Get z-values
  t_lam = t_lam.*n_plies
  z0 = sum(t_lam)/2.0
  z = zeros(n_lamina+1)-z0
  for k = 1:n_lamina
    z[k+1] = z[k+1]+sum(t_lam[1:k])
  end

  sc = S\[Nx,Ny,Nxy,Mx,My,Mxy]

  # Find strain at top and bottom of each lamina
  epsx = zeros(n_lamina+1)
  epsy = zeros(n_lamina+1)
  gammaxy = zeros(n_lamina+1)
  for k = 1:n_lamina+1
    epsx[k] = sc[1]+z[k]*sc[4]
    epsy[k] = sc[2]+z[k]*sc[5]
    gammaxy[k] = sc[3]+z[k]*sc[6]
  end

  # Rotate to align with material principal axes
  eps1U = zeros(n_lamina)
  eps1L = zeros(n_lamina)
  eps2U = zeros(n_lamina)
  eps2L = zeros(n_lamina)
  gamma12U = zeros(n_lamina)
  gamma12L = zeros(n_lamina)
  for k = 1:n_lamina
    theta = -tht_lam[k]
    strain = rotstrain(theta,epsx[k],epsy[k],gammaxy[k])
    eps1L[k] = strain[1]
    eps2L[k] = strain[2]
    eps12L[k] = strain[3]
    strain = rotstrain(theta,epsx[k+1],epsy[k+1],gammaxy[k+1])
    eps1U[k] = strain[1]
    eps2U[k] = strain[2]
    eps12U[k] = strain[3]
  end

return epsxU,epsxL,epsyU,epsyL,gammaxyU,gammaxyL
end

"""
    plystress(S,Nx,Ny,Nxy,Mx,My,Mxy,e1,e2,g12,anu12,n_lamina,n_plies,t_lam,
              tht_lam, mat_lam)
Calculates stress in each ply aligned with principal material direction
Returns sigma1U,sigma1L,sigma2U,sigma2L,tau12U,tau12L
"""
function plystress(S,Nx,Ny,Nxy,Mx,My,Mxy,e1::Array{Float64,1}, e2::Array{Float64,1},
  g12::Array{Float64,1}, anu12::Array{Float64,1},
  n_lamina::Int64,n_plies::Array{Int64,1}, t_lam::Array{Float64,1},
  tht_lam::Array{Float64,1}, mat_lam::Array{Int64,1})

  # get ply strains aligned with material principal axes
  eps1U,eps1L,eps2U,eps2L,gamma12U,gamma12L = plystrain(S,Nx,Ny,Nxy,Mx,My,Mxy,
                                              n_lamina,n_plies,t_lam,tht_lam)

  # material properties
  anud = 1.0 - anu12.*anu12.*e2./e1
  q11 = e1 ./ anud
  q22 = e2 ./ anud
  q12 = anu12.*e2 ./ anud
  q66 = g12

  # Convert strains to stresses
  sigma1U = zeros(n_lamina)
  sigma1L = zeros(n_lamina)
  sigma2U = zeros(n_lamina)
  sigma2L = zeros(n_lamina)
  tau12U = zeros(n_lamina)
  tau12L = zeros(n_lamina)
  for k = 1:n_lamina
    matid = mat_lam[k]
    sigma1U[k] = q11[matid]*eps1U[k] + q12[matid]*eps2U[k]
    sigma2U[k] = q12[matid]*eps1U[k] + q22[matid]*eps2U[k]
    tau12U[k] = q66[matid]*gamma12U[k]
    sigma1L[k] = q11[matid]*eps1L[k]+q12[matid]*eps2L[k]
    sigma2L[k] = q12[matid]*eps1L[k]+q22[matid]*eps2L[k]
    tau12L[k] = q66[matid]*gamma12L[k]
  end

  return sigma1U,sigma1L,sigma2U,sigma2L,tau12U,tau12L
end

"""
    rotmatprop(theta,qxx,qyy,qxy,qss)

Rotates material stiffness matrix theta degrees
"""
function rotmatprop(theta,qxx,qyy,qxy,qss)
  m = cosd(theta)
  n = sind(theta)
  qt = zeros(3,3)
  qt[1,1] = m^4.0*qxx + n^4.0*qyy + 2.0*m^2.0*n^2.0*qxy + 4.0*m^2.0*n^2.0*qss
  qt[1,2] = m^2.0*n^2.0*qxx + m^2.0*n^2.0*qyy+(m^4.0+n^4.0)*qxy - m^2.0*n^2.0*qss
  qt[1,3] = m^3.0*n*qxx - m*n^3.0*qyy + (m*n^3.0-m^3.0*n)*qxy + 2.0*(m*n^3.0-m^3.0*n)*qss
  qt[2,1] = qt[1,2]
  qt[2,2] = n^4.0*qxx + m^4.0*qyy + 2.0*m^2.0*n^2.0*qxy + 4.0*m^2.0*n^2.0*qss
  qt[2,3] = m*n^3.0*qxx - m^3.0*n*qyy + (m^3.0*n-m*n^3.0)*qxy + 2.0*(m^3.0*n-m*n^3.0)*qss
  qt[3,1] = qt[1,3]
  qt[3,2] = qt[2,3]
  qt[3,3] = m^2.0*n^2.0*qxx + m^2.0*n^2.0*qyy-2.0*m^2.0*n^2.0*qxy + (m^2.0-n^2.0)^2.0*qss
  return qt
end

"""
    rotstress(theta,sigmax,sigmay,tauxy)

Rotates stress theta degrees
"""
function rotstress(theta,sigmax,sigmay,tauxy)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0 st^2.0 2.0*st*ct; st^2.0 ct^2.0 -2.0*st*ct; -st*ct st*ct ct^2.0-st^2.0]
  return T*[sigmax,sigmay,tauxy]
end

"""
    rotstrain(theta,epsx,epsy,gammaxy)

Rotates strain theta degrees.
"""
function rotstrain(theta,epsx,epsy,gammaxy)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0 st^2.0 st*ct; st^2.0 ct^2.0 -st*ct; -2.0*st*ct 2.0*st*ct ct^2.0-st^2.0]
  return T*[epsx,epsy,gammaxy]
end

end #module
