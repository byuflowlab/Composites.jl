"""
    localbuckling(A::Array{<:Real,2}, D::Array{<:Real,2}, b::Real)
Calculates local buckling load and strains. `b` is the panel width.

assumes: balanced, symmetric, large aspect ratio, simply supported,
uniaxial compression, flat rectangular plate. buckling strain is output as a
positive value

see chapter on Structural Component Design Techniques from Alastair Johnson
section 6.2: Design of composite panels
"""
function localbuckling(A::AbstractArray{<:Real,2}, D::AbstractArray{<:Real,2},
    b::Real)
    bucklingload = 2.0*(pi/b)^2.0*(sqrt(D[1,1]*D[2,2])+D[1,2]+2*D[3,3])
    bucklingstrain = bucklingload/(A[1,1]-A[1,2]^2.0/A[2,2])
    return bucklingload, bucklingstrain
end #buckling
