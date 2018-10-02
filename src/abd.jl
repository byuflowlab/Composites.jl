"""
    `getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1})`

    `getz(lam::laminate)`

Returns a laminate's z-coordinates (coordinates of top and bottom of laminas)
given the thickness of plies in each lamina and the number of plies in each
lamina.
"""
function getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1})
    tlam = tply.*nply
    return pushfirst!(cumsum(tlam), 0.0).-sum(tlam)/2.0
end

getz(lam::laminate) = getz(lam.tply, lam.nply)

"""
    `getQ(e1::Real, e2::Real, g12::Real, nu12::Real, theta::Real=0.0)`

    `getQ(mat::material, theta::Real)`

    `getQ(mat::material, lam::laminate)`

Returns Q matrix rotated `theta` degrees
"""
function getQ(e1::Real, e2::Real, g12::Real, nu12::Real, theta::Real=0.0)
    nud = 1.0-e2/e1*nu12^2.0
    q11 = e1/nud
    q22 = e2/nud
    q12 = nu12*e2/nud
    q66 = g12
    return rotQ(q11, q12, q22, q66, theta)
end

getQ(mat::material, theta::Real=0.0) = getQ(mat.e1, mat.e2, mat.g12, mat.nu12,
    theta)

getQ(mat::Array{material,1}, lam::laminate) = getQ.(mat, lam.theta)

"""
    `rotQ(q::AbstractArray{<:Real,2}, theta::Real)`

    `rotQ(qxx::Real, q12::Real, q22::Real, q66::Real, theta::Real)`

Rotates Q matrix by theta degrees
"""
function rotQ(q::AbstractArray{<:Real,2}, theta::Real)
    if theta != 0.0
        c = cosd(theta)
        s = sind(theta)
        tsigma = [c^2.0 s^2.0 2*c*s; s^2.0 c^2.0 -2*c*s; -c*s c*s c^2-s^2]
        teps = LinearAlgebra.Diagonal([1.0,1.0,2.0])*tsigma*LinearAlgebra.Diagonal([1.0,1.0,0.5])
        qbar = tsigma\q*teps
    else
        qbar = q
    end
    return qbar
end

function rotQ(q11::Real, q12::Real, q22::Real, q66::Real, theta::Real)
  q = [q11 q12 0.0; q12 q22 0.0; 0.0 0.0 q66]
  return rotQ(q, theta)
end

"""
    `getABD(matid::AbstractArray{<:Integer,1},
        nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
        theta::AbstractArray{<:Real,1},
        q::AbstractArray{<:AbstractArray{<:Real,2},1})`

    `getABD(lam::laminate, q::AbstractArray{<:AbstractArray{<:Real,2},1})`

Returns A, B, and D matrices

# Arguments:
* `matid::AbstractArray{<:Integer,1}`: material id of each lamina
* `nply::AbstractArray{<:Integer,1}`: number of plies in each lamina
* `tply::AbstractArray{<:Real,1}`: thickness of a ply (m) in each lamina
* `theta::AbstractArray{<:Real,1}`: orientation (deg) of each lamina
* `q::AbstractArray{<:AbstractArray{<:Real,2}}`: Stiffness matrix of each lamina
"""
function getABD(matid::AbstractArray{<:Integer,1},
    nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, q::AbstractArray{<:AbstractArray{<:Real,2},1})

    R = typeof(q[1][1])

    nlam = length(nply)

    z = getz(tply, nply)

    # Loop through layers filling in ABD matrix
    A = zeros(R, 3, 3)
    B = zeros(R, 3, 3)
    D = zeros(R, 3, 3)
    for k = 1:nlam
        # Rotate material stiffness properties by specifed angle
        imat = matid[k]
        qbar = rotQ(q[imat], theta[k])

        # Find lumped stiffness properties for the laminate
        for i = 1:3
            for j = 1:3
                A[i,j] = A[i,j]+qbar[i,j]*(z[k+1]-z[k])
                B[i,j] = B[i,j]+qbar[i,j]/2.0*(z[k+1]^2.0-z[k]^2.0)
                D[i,j] = D[i,j]+qbar[i,j]/3.0*(z[k+1]^3.0-z[k]^3.0)
            end
        end
    end

    return A, B, D
end
getABD(lam::laminate, q::AbstractArray{<:AbstractArray{<:Real,2},1}) =
    getABD(lam.matid, lam.nply, lam.tply, lam.theta, q)

"""
    compliance(A::Array{<:Real,2}, B::Array{<:Real,2}, D::Array{<:Real,2})
Returns alpha, beta, and delta
"""
function compliance(A::AbstractArray{<:Real,2}, B::AbstractArray{<:Real,2},
    D::AbstractArray{<:Real,2})
    invA = inv(A)
    delta = inv(D-B*invA*B)
    beta = -invA*B*delta
    alfa = invA-beta*B*invA
    return alfa, beta, delta
end # compliance
