"""
    getz(tply::Array{Real,1},nply::Array{Int64,1})
    getz(lam::laminate)
Returns a laminate's z-coordinates (coordinates of top and bottom of laminas)
given the thickness of plies in each lamina and the number of plies in each
lamina.
"""
function getz(tply::Array{Real,1},nply::Array{Int64,1})
  # Get z-values
  tlam = tply.*nply
  return unshift!(cumsum(tlam),0.0)-sum(tlam)/2.0
end
getz(lam::laminate) = getz(lam.tply,lam.nply)


"""
    getQ(e1::Real,e2::Real,g12::Real,nu12::Real,theta::Real=0.0)
    getQ(mat::material,theta::Real)
Returns Q matrix rotated by angle theta
"""
function getQ(e1::Real,e2::Real,g12::Real,nu12::Real,theta::Real=0.0)
  nud = 1.0-e2/e1*nu12^2.0
  q11 = e1/nud
  q22 = e2/nud
  q12 = nu12*e2/nud
  q66 = g12
  return rotQ(q11,q12,q22,q66,theta)
end
getQ(mat::material,theta::Real=0.0) = getQ(mat.e1,mat.e2,mat.g12,mat.nu12,theta)

"""
    getQ(e1::Array{Real,1},e2::Array{Real,1},g12::Array{Real,1},
      nu12::Array{Real,1},theta::Array{Real,1}=zeros(Real,length(e1)))
    getQ(mat::Array{material,1},theta::Array{Real,1}=zeros(Real,length(mat)))
    getQ(mat::Array{material,1},lam::laminate)
Returns multiple Q matrices. Third dimension corresponds to the ply number.
"""
function getQ(e1::Array{Real,1},e2::Array{Real,1},g12::Array{Real,1},
  nu12::Array{Real,1},theta::Array{Real,1}=zeros(Real,length(e1)))
  nply = length(e1)
  qbar = zeros(Real,3,3,nply)
  for i = 1:nply
    qbar[:,:,i] = getQ(e1[i],e2[i],g12[i],nu12[i],theta[i])
  end
  return qbar
end

function getQ(mat::Array{material,1},theta::Array{Real,1}=zeros(Real,length(mat)))
  nply = length(mat)
  qbar = zeros(Real,3,3,nply)
  for i = 1:nply
    qbar[:,:,i] = getQ(mat[i].e1,mat[i].e2,mat[i].g12,mat[i].nu12,theta[i])
  end
  return qbar
end
getQ(mat::Array{material,1},lam::laminate) = getQ(mat,lam.theta)

"""
    rotQ(qxx::Real,q12::Real,q22::Real,q66::Real,theta::Real)
Rotates Q matrix by theta degrees
"""
function rotQ(q11::Real,q12::Real,q22::Real,q66::Real,theta::Real)
  q = [q11 q12 0.0;q12 q22 0.0; 0.0 0.0 q66]
  return rotQ(q,theta)
end

"""
    rotQ(q::Array{Real,2},theta::Real)
Rotates Q matrix by theta degrees
"""
function rotQ(q::Array{Real,2},theta::Real)
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
    getABD(nply::Array{Int64,1},tply::Array{Real,1},matid::Array{Int64,1},
    theta::Array{Real,1},e1::Array{Real,1},e2::Array{Real,1},
    g12::Array{Real,1},nu12::Array{Real,1})
    getABD(lam::laminate,q::Array{Real,3}) = getABD(lam.matid,lam.nply,lam.tply,lam.theta,q)
Returns A, B, and D matrices
# Arguments:
* `matid::Array{Int64,1}`: material id of each lamina
* `nply::Array{Int64,1}`: number of plies in each lamina
* `tply::Array{Real,1}`: thickness of a ply (m) in each lamina
* `theta::Array{Real,1}`: orientation (deg) of each lamina
* `e1::Array{Real,1}`: E1 corresponding to material id
* `e2::Array{Real,1}`: E2 corresponding to material id
* `g12::Array{Real,1}`: G12 corresponding to material id
* `nu12::Array{Real,1}`: nu12 corresponding to material id
"""
function getABD(matid::Array{Int64,1},nply::Array{Int64,1},tply::Array{Real,1},
  theta::Array{Real,1},q::Array{Real,3})

  nlam = length(nply)

  z = getz(tply,nply)

  # Loop through layers filling in ABD matrix
  A = zeros(3,3)
  B = zeros(3,3)
  D = zeros(3,3)
  for k = 1:nlam
    # Rotate material stiffness properties by specifed angle
    imat = matid[k]
    qbar = rotQ(q[:,:,imat],theta[k])

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
getABD(lam::laminate,q::Array{Real,3}) = getABD(lam.matid,lam.nply,lam.tply,lam.theta,q)

"""
    compliance(A::Array{Real,2},B::Array{Real,2},D::Array{Real,2})
Returns alpha, beta, and delta
"""
function compliance(A::Array{Real,2},B::Array{Real,2},D::Array{Real,2})
  invA = inv(A)
  delta = inv(D-B*invA*B)
  beta = -invA*B*delta
  alfa = invA-beta*B*invA
  return alfa,beta,delta
end # compliance
