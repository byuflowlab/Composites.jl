"""
    `getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1})`

    `getz(lam::Laminate)`

Returns a laminate's z-coordinates (coordinates of top and bottom of laminas)
given the thickness of plies in each lamina and the number of plies in each
lamina.
"""
function getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1})
    ttotal = 0.0
    tlam = zeros(eltype(tply), length(tply)+1)
    for i = 1:length(tply)
        ttotal += tply[i]*nply[i]
        tlam[i+1] = ttotal
    end
    tlam .-= ttotal/2.0
    return tlam
end
getz(lam::Laminate) = getz(lam.tply, lam.nply)

"""
    `getQ(e1::Real, e2::Real, g12::Real, nu12::Real, theta::Real=0.0)`

    `getQ(mat::Material, theta::Real)`

    `getQ(mat::Material, lam::Laminate)`

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
getQ(mat::Material, theta::Real=0.0) = getQ(mat.e1, mat.e2, mat.g12, mat.nu12, theta)
getQ(mat::Array{Material,1}, lam::Laminate) = getQ.(mat, lam.theta)

"""
    `rotQ(q::AbstractArray{<:Real,2}, theta::Real)`

    `rotQ(qxx::Real, q12::Real, q22::Real, q66::Real, theta::Real)`

Rotates Q matrix by theta degrees
"""
function rotQ!(q::AbstractArray{<:Real,2}, theta::Real)

    if theta != zero(eltype(theta))
        c = cosd(theta)
        s = sind(theta)
        c2 = c^2
        s2 = s^2
        cs = c*s

        # qbar = inv(tsigma)*q*teps

        # tsigma = [c^2.0 s^2.0 2*c*s; s^2.0 c^2.0 -2*c*s; -c*s c*s c^2-s^2]

        inv_tsig11 = c2
        inv_tsig12 = s2
        inv_tsig13 = -2*cs
        inv_tsig21 = s^2
        inv_tsig22 = c^2
        inv_tsig23 = 2*cs
        inv_tsig31 = cs
        inv_tsig32 = -cs
        inv_tsig33 = 2*c2-1

        # teps = LinearAlgebra.Diagonal([1.0,1.0,2.0])*tsigma*LinearAlgebra.Diagonal([1.0,1.0,0.5])

        teps11 = c2
        teps12 = s2
        teps13 = cs
        teps21 = s2
        teps22 = c2
        teps23 = -cs
        teps31 = -2*cs
        teps32 = 2*cs
        teps33 = c2-s2

        # inv(tsigma)*q

        inv_tsig_q11 = inv_tsig11*q[1,1] + inv_tsig12*q[2,1] + inv_tsig13*q[3,1]
        inv_tsig_q12 = inv_tsig11*q[1,2] + inv_tsig12*q[2,2] + inv_tsig13*q[3,2]
        inv_tsig_q13 = inv_tsig11*q[1,3] + inv_tsig12*q[2,3] + inv_tsig13*q[3,3]
        inv_tsig_q21 = inv_tsig21*q[1,1] + inv_tsig22*q[2,1] + inv_tsig23*q[3,1]
        inv_tsig_q22 = inv_tsig21*q[1,2] + inv_tsig22*q[2,2] + inv_tsig23*q[3,2]
        inv_tsig_q23 = inv_tsig21*q[1,3] + inv_tsig22*q[2,3] + inv_tsig23*q[3,3]
        inv_tsig_q31 = inv_tsig31*q[1,1] + inv_tsig32*q[2,1] + inv_tsig33*q[3,1]
        inv_tsig_q32 = inv_tsig31*q[1,2] + inv_tsig32*q[2,2] + inv_tsig33*q[3,2]
        inv_tsig_q33 = inv_tsig31*q[1,3] + inv_tsig32*q[2,3] + inv_tsig33*q[3,3]

        # inv(tsigma)*q*teps

        q[1,1] = inv_tsig_q11*teps11 + inv_tsig_q12*teps21 + inv_tsig_q13*teps31
        q[1,2] = inv_tsig_q11*teps12 + inv_tsig_q12*teps22 + inv_tsig_q13*teps32
        q[1,3] = inv_tsig_q11*teps13 + inv_tsig_q12*teps23 + inv_tsig_q13*teps33
        q[2,1] = inv_tsig_q21*teps11 + inv_tsig_q22*teps21 + inv_tsig_q23*teps31
        q[2,2] = inv_tsig_q21*teps12 + inv_tsig_q22*teps22 + inv_tsig_q23*teps32
        q[2,3] = inv_tsig_q21*teps13 + inv_tsig_q22*teps23 + inv_tsig_q23*teps33
        q[3,1] = inv_tsig_q31*teps11 + inv_tsig_q32*teps21 + inv_tsig_q33*teps31
        q[3,2] = inv_tsig_q31*teps12 + inv_tsig_q32*teps22 + inv_tsig_q33*teps32
        q[3,3] = inv_tsig_q31*teps13 + inv_tsig_q32*teps23 + inv_tsig_q33*teps33

    end

    return q
end

function rotQ(q::AbstractArray{<:Real,2}, theta::Real)
    qbar = copy(q)
    return rotQ!(qbar, theta)
end

function rotQ(q11::Real, q12::Real, q22::Real, q66::Real, theta::Real)
  q0 = zero(eltype(q11))
  q = [q11 q12 q0; q12 q22 q0; q0 q0 q66]
  return rotQ!(q, theta)
end

"""
    `getABD(matid::AbstractArray{<:Integer,1},
        nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
        theta::AbstractArray{<:Real,1},
        q::AbstractArray{<:AbstractArray{<:Real,2},1})`

    `getABD(lam::Laminate, q::AbstractArray{<:AbstractArray{<:Real,2},1})`

Returns A, B, and D matrices

# Arguments:
* `matid::AbstractArray{<:Integer,1}`: material id of each lamina
* `nply::AbstractArray{<:Integer,1}`: number of plies in each lamina
* `tply::AbstractArray{<:Real,1}`: thickness of a ply (m) in each lamina
* `theta::AbstractArray{<:Real,1}`: orientation (deg) of each lamina
* `q::AbstractArray{<:AbstractArray{<:Real,2}}`: Stiffness matrix of each lamina
"""
function getABD!(A::AbstractArray{<:Real,2}, B::AbstractArray{<:Real,2}, D::AbstractArray{<:Real,2},
    matid::AbstractArray{<:Integer,1}, nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, q::AbstractArray{<:AbstractArray{<:Real,2},1})

    typeof(A)

    R = eltype(tply)

    nlam = length(nply)

    z = getz(tply, nply)

    qbar = zeros(R, 3, 3)

    # Loop through layers filling in ABD matrix
    for k = 1:nlam
        # Rotate material stiffness properties by specifed angle
        imat = matid[k]
        qbar .= q[imat]
        rotQ!(qbar, theta[k])

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

function getABD(matid::AbstractArray{<:Integer,1}, nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, q::AbstractArray{<:AbstractArray{<:Real,2},1})

    R = eltype(tply)

    A = zeros(R, 3, 3)
    B = zeros(R, 3, 3)
    D = zeros(R, 3, 3)

    return getABD!(A, B, D, matid, nply, tply, theta, q)
end

getABD(lam::Laminate, q::AbstractArray{<:AbstractArray{<:Real,2},1}) =
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
