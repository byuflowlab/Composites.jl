# This needs to be rewritten
# Modified 10 October 2017 TGM
"""
    localbuckling(A::Array{<:Real,2},D::Array{<:Real,2},b::Real)
Calculates local buckling load and strains

assumes: balanced, symmetric,large aspect ratio, simply supported,
uniaxial compression, flat rectangular plate

see chapter on Structural Component Design Techniques from Alastair Johnson
section 6.2: Design of composite panels
"""
function localbuckling(A::Array{<:Real,2},D::Array{<:Real,2},b::Real)
  bucklingload = 2.0*(pi/b)^2.0*(sqrt(D[1,1]*D[2,2])+D[1,2]+2*D[3,3])
  bucklingstrain = -bucklingload/(A[1,1]-A[1,2]^2.0/A[2,2])
  return bucklingload,bucklingstrain
end #buckling

function localbuckling(A::Array{Array{Real,2},1},D::Array{Array{Real,2},1},b::Array{<:Real,1})
  bucklingload = zeros(Real,length(A))
  bucklingstrain = zeros(Real,length(A))
  for i = 1:length(A)
    bucklingload[i],bucklingstrain[i] = localbuckling(A[i],D[i],b[i])
  end
  return bucklingload,bucklingstrain
end
