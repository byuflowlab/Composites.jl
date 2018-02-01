"""
    getQ(e1::Float64,e2::Float64,g12::Float64,nu12::Float64,theta::Float64=0.0)
Returns Q matrix rotated by angle theta
"""
function getQ(e1::Float64,e2::Float64,g12::Float64,nu12::Float64,theta::Float64=0.0)
  nud = 1.0-e2/e1*nu12^2.0
  q11 = e1/nud
  q22 = e2/nud
  q12 = nu12*e2/nud
  q66 = g12
  return rotQ(q11,q12,q22,q66,theta)
end

"""
    plystiffness(e1::Array{Float64,1},e2::Array{Float64,1},g12::Array{Float64,1},nu12::Array{Float64,1},theta::Array{Float64,1}=0.0)
Returns multiple Q matrices. Third dimension corresponds to the ply number.
"""
function getQ(e1::Array{Float64,1},e2::Array{Float64,1},g12::Array{Float64,1},nu12::Array{Float64,1},theta::Array{Float64,1}=zeros(Float64,length(e1)))
  # Implicitly assigned inputs
  if !(length(e1) == length(e2) == length(g12) == length(nu12) == length(theta))
    error("lengths of specified material properties do not match")
  end
  nply = length(e1)
  nud = 1.0-e2./e1.*nu12.^2.0
  q11 = e1./nud
  q22 = e2./nud
  q12 = nu12.*e2./nud
  q66 = g12
  qbar = zeros(Float64,3,3,nply)
  for i = 1:nply
    qbar[:,:,i] = rotQ(q11[i],q12[i],q22[i],q66[i],theta[i])
  end
  return qbar
end

"""
    rotQ(qxx::Float64,q12::Float64,q22::Float64,q66::Float64,theta::Float64)
Rotates Q matrix by theta degrees
"""
function rotQ(q11::Float64,q12::Float64,q22::Float64,q66::Float64,theta::Float64)
  q = [q11 q12 0.0;q12 q22 0.0; 0.0 0.0 q66]
  return rotQ(q,theta)
end

"""
    rotQ(q::Array{Float64,2},theta::Float64)
Rotates Q matrix by theta degrees
"""
function rotQ(q::Array{Float64,2},theta::Float64)
  if theta != 0.0
    c = cosd(theta)
    s = sind(theta)
    tsigma = [c^2.0 s^2.0 2*c*s;s^2.0 c^2.0 -2*c*s;-c*s c*s c^2-s^2]
    teps = Diagonal([1.0,1.0,2.0])*tsigma*Diagonal([1.0,1.0,0.5])
    qbar = tsigma\q*teps
  else
    qbar = q
  end
  return qbar
end

"""
    getz(tply::Array{Float64,1},nply::Array{Int64,1})
Returns a laminate's z-coordinates (coordinates of top and bottom of laminas)
given the thickness of plies in each lamina and the number of plies in each
lamina.
"""
function getz(tply::Array{Float64,1},nply::Array{Int64,1})
  # Get z-values
  tlam = tply.*nply
  return unshift!(cumsum(tlam),0.0)-sum(tlam)/2.0
end

"""
    getABD(nply::Array{Int64,1},tply::Array{Float64,1},matlam::Array{Int64,1},
    thetalam::Array{Float64,1},e1::Array{Float64,1},e2::Array{Float64,1},
    g12::Array{Float64,1},nu12::Array{Float64,1})
Returns A, B, and D matrices
# Arguments:
* `nply::Array{Int64,1}`: number of plies in each lamina
* `tply::Array{Float64,1}`: thickness of a ply (m) in each lamina
* `matlam::Array{Int64,1}`: material id of each lamina
* `thetalam::Array{Float64,1}`: orientation (deg) of each lamina
* `e1::Array{Float64,1}`: E1 corresponding to material id
* `e2::Array{Float64,1}`: E2 corresponding to material id
* `g12::Array{Float64,1}`: G12 corresponding to material id
* `nu12::Array{Float64,1}`: nu12 corresponding to material id
"""
function getABD(nply::Array{Int64,1},tply::Array{Float64,1},matlam::Array{Int64,1},
  thetalam::Array{Float64,1},e1::Array{Float64,1},e2::Array{Float64,1},
  g12::Array{Float64,1},nu12::Array{Float64,1})

  nlam = length(nply)

  z = getz(tply,nply)

  # Get Q matrices
  q = getQ(e1,e2,g12,nu12)

  # Loop through layers filling in ABD matrix
  A = zeros(3,3)
  B = zeros(3,3)
  D = zeros(3,3)
  for k = 1:nlam
    # Rotate material stiffness properties by specifed angle
    imat = matlam[k]
    qbar = rotQ(q[:,:,imat],thetalam[k])

    # Find lumped stiffness properties for the laminate
    for i = 1:3
      for j = 1:3
        A[i,j] = A[i,j]+qbar[i,j]*(z[k+1]-z[k])
        B[i,j] = B[i,j]+qbar[i,j]/2.0*(z[k+1]^2.0-z[k]^2.0)
        D[i,j] = D[i,j]+qbar[i,j]/3.0*(z[k+1]^3.0-z[k]^3.0)
      end
    end
  end

  return A,B,D
end #stiffness

"""
    compliance(A::Array{Float64,2},B::Array{Float64,2},D::Array{Float64,2})
Returns alpha, beta, and delta
"""
function compliance(A::Array{Float64,2},B::Array{Float64,2},D::Array{Float64,2})
  invA = inv(A)
  delta = inv(D-B*invA*B)
  beta = -invA*B*delta
  alfa = invA-beta*B*invA
  return alfa,beta,delta
end # compliance
