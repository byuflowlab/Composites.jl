struct material
  e1::T where T<:Real #E1
  e2::T where T<:Real #E2
  g12::T where T<:Real #G12
  nu12::T where T<:Real #nu12
  rho::T where T<:Real #density
  xt::T where T<:Real #longitudinal tensile ultimate strength (Pa)
  xc::T where T<:Real #longitudinal compressive ultimate strength (Pa)
  yt::T where T<:Real #transverse tensile ultimate strength (Pa)
  yc::T where T<:Real #transverse compressive ultimate strength (Pa)
  s::T where T<:Real #shear ultimate strength (Pa)
end

function material(e1::Array{Real,1},e2::Array{Real,1},g12::Array{Real,1},
  nu12::Array{Real,1},rho::Array{Real,1},xt::Array{Real,1},
  xc::Array{Real,1},yt::Array{Real,1},yc::Array{Real,1},
  s::Array{Real,1})

  nmat = length(e1)
  matprops = Array{material,1}(nmat)
  for i = 1:nmat
    matprops[i] = material(e1[i],e2[i],g12[i],nu12[i],rho[i],xt[i],
                           xc[i],yt[i],yc[i],s[i])
  end
  return matprops
end

struct laminate
  matid::Array{Int64,1}
  nply::Array{Int64,1}
  tply::Array{T,1} where T<:Real
  theta::Array{T,1} where T<:Real
  function laminate(matid::Array{Int64,1},nply::Array{Int64,1},tply::Array{Real,1},theta::Array{Real,1})
    if !(length(matid) == length(nply) == length(tply) == length(theta))
      error("All arrays must have same length")
    end
    return new(matid,nply,tply,theta)
  end
end
