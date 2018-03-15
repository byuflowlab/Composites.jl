struct material{R<:Real}
  e1::R #E1
  e2::R #E2
  g12::R #G12
  nu12::R #nu12
  rho::R #density
  xt::R #longitudinal tensile ultimate strength (Pa)
  xc::R #longitudinal compressive ultimate strength (Pa)
  yt::R #transverse tensile ultimate strength (Pa)
  yc::R #transverse compressive ultimate strength (Pa)
  s::R #shear ultimate strength (Pa)
end

function material(e1::Array{<:Real,1},e2::Array{<:Real,1},g12::Array{<:Real,1},
  nu12::Array{<:Real,1},rho::Array{<:Real,1},xt::Array{<:Real,1},
  xc::Array{<:Real,1},yt::Array{<:Real,1},yc::Array{<:Real,1},
  s::Array{<:Real,1})

  nmat = length(e1)
  matprops = Array{material,1}(nmat)
  for i = 1:nmat
    matprops[i] = material(e1[i],e2[i],g12[i],nu12[i],rho[i],xt[i],
                           xc[i],yt[i],yc[i],s[i])
  end
  return matprops
end

struct laminate{I<:Integer,R<:Real}
  matid::Array{I,1}
  nply::Array{I,1}
  tply::Array{R,1}
  theta::Array{R,1}
  function laminate{I<:Integer,R<:Real}(matid::Array{<:Integer,1},nply::Array{<:Integer,1},tply::Array{<:Real,1},theta::Array{<:Real,1})
    if !(length(matid) == length(nply) == length(tply) == length(theta))
      error("All arrays must have same length")
    end
    return new(matid,nply,tply,theta)
  end
end
