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
  bucklingstrain = -bucklingload/(A[1,1]-A[1,2]^2.0/A[2,2])
  return bucklingload,bucklingstrain
end #buckling

function localbuckling(A::Array{Array{Float64,2},1},D::Array{Array{Float64,2},1},b::Array{Float64,1})
  bucklingload = zeros(Float64,length(A))
  bucklingstrain = zeros(Float64,length(A))
  for i = 1:length(A)
    bucklingload[i],bucklingstrain[i] = localbuckling(A[i],D[i],b[i])
  end
  return bucklingload,bucklingstrain
end
