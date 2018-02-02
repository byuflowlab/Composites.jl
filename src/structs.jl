struct material
  e1::Float64 #E1
  e2::Float64 #E2
  g12::Float64 #G12
  nu12::Float64 #nu12
  rho::Float64 #density
  xt::Float64 #longitudinal tensile ultimate strength (Pa)
  xc::Float64 #longitudinal compressive ultimate strength (Pa)
  yt::Float64 #transverse tensile ultimate strength (Pa)
  yc::Float64 #transverse compressive ultimate strength (Pa)
  s::Float64 #shear ultimate strength (Pa)
end

function material(e1::Array{Float64,1},e2::Array{Float64,1},g12::Array{Float64,1},
  nu12::Array{Float64,1},rho::Array{Float64,1},xt::Array{Float64,1},
  xc::Array{Float64,1},yt::Array{Float64,1},yc::Array{Float64,1},
  s::Array{Float64,1})

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
  tply::Array{Float64,1}
  theta::Array{Float64,1}
  function laminate(matid::Array{Int64,1},nply::Array{Int64,1},tply::Array{Float64,1},theta::Array{Float64,1})
    if !(length(matid) == length(nply) == length(tply) == length(theta))
      error("All arrays must have same length")
    end
    return new(matid,nply,tply,theta)
  end
end
