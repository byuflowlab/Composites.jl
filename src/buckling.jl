# This needs to be rewritten
# Modified 10 October 2017 TGM
"""
    localbuckling(A::Array{Float64,2},D::Array{Float64,2},b::Float64)
Calculates local buckling load and strains

assumes: balanced, symmetric,large aspect ratio, simply supported,
uniaxial compression, flat rectangular plate

see chapter on Structural Component Design Techniques from Alastair Johnson
section 6.2: Design of composite panels
"""
function localbuckling(A::Array{Float64,2},D::Array{Float64,2},b::Float64)
  bucklingload = 2.0*(pi/b)^2.0*(sqrt(D[1,1]*D[2,2])+D[1,2]+2*D[3,3])
  bucklingstrain = -Ncr1/(A[1,1]-A[1,2]^2.0/A[2,2])
  return bucklingload,bucklingstrain
end #buckling
