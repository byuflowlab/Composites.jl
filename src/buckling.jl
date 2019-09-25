"""
    plate_uniaxial_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, b::Real)
Calculates local buckling load and strains. `b` is the panel width.

assumes: balanced, symmetric, large aspect ratio, simply supported,
uniaxial compression, flat rectangular plate. buckling strain is output as a
positive value

see chapter on Structural Component Design Techniques from Alastair Johnson
section 6.2: Design of composite panels

probably assumes wave numbers are continuous in order to eliminate them
"""
function plate_uniaxial_local_buckling(A::AbstractArray{<:Real,2}, D::AbstractArray{<:Real,2},
    b::Real)
    bucklingload = 2.0*(pi/b)^2.0*(sqrt(D[1,1]*D[2,2])+D[1,2]+2*D[3,3])
    # bucklingstrain = bucklingload/(A[1,1]-A[1,2]^2.0/A[2,2])
    return bucklingload
end #buckling

# Lekhnitskii (see LOCAL BUCKLING OF COMPOSITE BEAMS G. Tarján, A. Sapkás and L.P. Kollár)
function plate_bending_local_buckling(D::AbstractArray{<:Real,2}, b::Real)
    bucklingload = (pi/b)^2.0*(13.4*sqrt(D[1,1]*D[2,2])+10.4*(D[1,2]+2*D[3,3]))
    return bucklingload
end #buckling

"""
    plate_shear_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, b::Real)
Calculates local buckling load and strains. `b` is the panel width.

assumes: balanced, symmetric, large aspect ratio, simply supported,
uniaxial compression, flat rectangular plate

see chapter on Structural Component Design Techniques from Alastair Johnson
section 6.2: Design of composite panels

probably assumes wave numbers are continuous in order to eliminate them
"""
function plate_shear_local_buckling(A::AbstractArray{<:Real,2}, D::AbstractArray{<:Real,2},
    b::Real)
    D1 = D[1,1]
    D2 = D[2,2]
    D3 = D[1,2] + 2*D[3,3]
    c = (D1*D2)^(1/2)/D3
    #TODO: replace if statement here with smooth version
    if c >= 1.0
        fc = c^(1/2)*(0.62+0.38/c)
    else
        fc = 0.89 + 0.04*c + 0.07*c^2
    end
    K2 = (D2*D3)^(1/2)*fc/D1
    bucklingload = 52*K2*D1/b^2
    # bucklingstrain = -bucklingload/A[3,3]
    return bucklingload
end #buckling

"""
    cylinder_uniaxial_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, R::Real, h::Real, correlation::Bool=true)
Calculates local buckling load. `R` is the cylinder radius. `h` is the total laminate thickness.

assumes: infinitely long cylinder, simply supported, uniaxial compression,
wave numbers are real-valued rather than integers (this allows them to be
eliminated as variables)

see An Approximate Solution and Master Curves for Buckling of Symmetrically
Laminated Composite Cylinders by Michael P. Nemeth

also see Simple Formulas and Results for Buckling Resistance and Stiffness Design of
Compression-Loaded Laminated-Composite Cylinders by Michael P. Nemeth

The correlation factor from SP-8001 is also applied, but can be disabled by
setting `correlation=false`.

"""
function cylinder_uniaxial_local_buckling(A, D, R, h, correlation=true)
    nu_m = A[1,2]/sqrt(A[1,1]*A[2,2])
    rho = R/h*sqrt(1-nu_m^2)*(A[1,1]*A[2,2]*h^4/(144*D[1,1]*D[2,2]))^(1/4)
    beta = (D[1,2] + 2*D[3,3])/sqrt(D[1,1]*D[2,2])
    I = (A[1,1]*D[2,2]/(D[1,1]*A[2,2]))^(1/2)
    pcr_tilde = 4*sqrt(3*I)*rho/(pi^2)
    bucklingload = pcr_tilde*pi^2*sqrt(D[1,1]*D[2,2])/(R^2)

    if correlation
        psi = 1/29.8*(R/(D[1,1]*D[2,2]/(A[1,1]*A[2,2]))^(1/4))^(1/2)
        gamma = 1-0.901*(1-exp(-psi))
        bucklingload = gamma*bucklingload
    end

    return bucklingload
end

"""
    cylinder_bending_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, R::Real, h::Real)
Calculates local buckling load. `R` is the cylinder radius. `h` is the total laminate thickness.

Follows approach suggested in SP-8001, which uses the theoretical axial compression cylinder
buckling equation with a correction factor specific for bending cases, except
assumes an infinitely long cylinder and that wave numbers are real-valued rather
than integers (in order to elimate them as variables)

see An Approximate Solution and Master Curves for Buckling of Symmetrically
Laminated Composite Cylinders by Michael P. Nemeth for details on the equation
for the buckling of an infinitely long orthotropic cylinder

"""
function cylinder_bending_local_buckling(A, D, R, h)

    bucklingload = cylinder_bending_local_buckling(A, D, R, h, correlation=false)

    # apply correlation
    psi = 1/29.8*(R/(D[1,1]*D[2,2]/(A[1,1]*A[2,2]))^(1/4))^(1/2)
    gamma = 1-0.731*(1-exp(-psi))
    bucklingload = gamma*bucklingload

    return bucklingload
end
