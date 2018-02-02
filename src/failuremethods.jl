# Alternative single evaluation methods
maxstress(stress::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = maxstress(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
maxstress(sigma1::Float64,sigma2::Float64,tau12::Float64,mat::material) =
  maxstress(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Alternative laminate evaluation methods
function maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Array{Float64,1},xc::Array{Float64,1},yt::Array{Float64,1},
  yc::Array{Float64,1},s::Array{Float64,1},imat::Array{Float64,1})
  nlam = length(imat)
  matfail = zeros(Float64,6,nlam)
  sf = zeros(Float64,6,nlam)
  for i = 1:nlam
    matfail[:,i],sf[:,i] = maxstress(sigma1[i],sigma2[i],tau12[i],xt[imat[i]],
      xc[imat[i]],yt[imat[i]],yc[imat[i]],s[imat[i]])
  end
  return matfail,sf
end
maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Array{Float64,1},xc::Array{Float64,1},yt::Array{Float64,1},
  yc::Array{Float64,1},s::Array{Float64,1},lam::laminate) = maxstress(sigma1,
  sigma2,tau12,xt,xc,yt,yc,s,lam.imat)
maxstress(stress::Array{Float64,2},xt::Array{Float64,1},xc::Array{Float64,1},
  yt::Array{Float64,1},yc::Array{Float64,1},s::Array{Float64,1},
  imat::Array{Float64,1}) = maxstress(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,imat)
maxstress(stress::Array{Float64,2},xt::Array{Float64,1},xc::Array{Float64,1},
  yt::Array{Float64,1},yc::Array{Float64,1},s::Array{Float64,1},
  lam::laminate) = maxstress(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,lam.imat)
# Alternative laminate evaluation methods using mat::Array{material,1}
function maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},mat::material,imat::Array{Float64,1})
  nlam = length(imat)
  matfail = zeros(Float64,6,nlam)
  for i = 1:nlam
    matfail[:,i] = maxstress(sigma1[i],sigma2[i],tau12[i],mat[imat[i]].xt,
      mat[imat[i]].xc,mat[imat[i]].yt,mat[imat[i]].yc,mat[imat[i]].s)
  end
  return matfail
end
maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},tau12::Array{Float64,1},
  mat::material,lam::laminate) = maxstress(sigma1,sigma2,tau12,mat,lam.imat)
maxstress(stress::Array{Float64,2},mat::material,imat::Array{Float64,1}) =
  maxstress(stress[1,:],stress[2,:],stress[3,:],mat,imat)
maxstress(stress::Array{Float64,2},mat::material,lam::laminate) =
  maxstress(stress[1,:],stress[2,:],stress[3,:],mat,lam.imat)

# Alternative single evaluation methods
tsaiwu(stress::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = tsaiwu(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
tsaiwu(sigma1::Float64,sigma2::Float64,tau12::Float64,mat::material) =
  tsaiwu(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Alternative laminate evaluation methods
function tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Array{Float64,1},xc::Array{Float64,1},yt::Array{Float64,1},
  yc::Array{Float64,1},s::Array{Float64,1},imat::Array{Float64,1})
  nlam = length(imat)
  matfail = zeros(Float64,nlam)
  sf = zeros(Float64,2,nlam)
  for i = 1:nlam
    matfail[i],sf[:,i] = tsaiwu(sigma1[i],sigma2[i],tau12[i],xt[imat[i]],
      xc[imat[i]],yt[imat[i]],yc[imat[i]],s[imat[i]])
  end
  return matfail,sf1,sf2
end
tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Array{Float64,1},xc::Array{Float64,1},yt::Array{Float64,1},
  yc::Array{Float64,1},s::Array{Float64,1},lam::laminate) = tsaiwu(sigma1,
  sigma2,tau12,xt,xc,yt,yc,s,lam.imat)
tsaiwu(stress::Array{Float64,2},xt::Array{Float64,1},xc::Array{Float64,1},
  yt::Array{Float64,1},yc::Array{Float64,1},s::Array{Float64,1},
  imat::Array{Float64,1}) = tsaiwu(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,imat)
tsaiwu(stress::Array{Float64,2},xt::Array{Float64,1},xc::Array{Float64,1},
  yt::Array{Float64,1},yc::Array{Float64,1},s::Array{Float64,1},
  lam::laminate) = tsaiwu(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,lam.imat)
# Alternative laminate evaluation methods using mat::Array{material,1}
function tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},mat::material,imat::Array{Float64,1})
  nlam = length(imat)
  matfail = zeros(Float64,nlam)
  for i = 1:nlam
    matfail[i] = tsaiwu(sigma1[i],sigma2[i],tau12[i],mat[imat[i]].xt,
      mat[imat[i]].xc,mat[imat[i]].yt,mat[imat[i]].yc,mat[imat[i]].s)
  end
  return matfail
end
tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},tau12::Array{Float64,1},
  mat::material,lam::laminate) = tsaiwu(sigma1,sigma2,tau12,mat,lam.imat)
tsaiwu(stress::Array{Float64,2},mat::material,imat::Array{Float64,1}) =
  tsaiwu(stress[1,:],stress[2,:],stress[3,:],mat,imat)
tsaiwu(stress::Array{Float64,2},mat::material,lam::laminate) =
  tsaiwu(stress[1,:],stress[2,:],stress[3,:],mat,lam.imat)


# Alternative single evaluation methods
hashinrotem(stress::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = hashinrotem(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
hashinrotem(sigma1::Float64,sigma2::Float64,tau12::Float64,mat::material) =
  hashinrotem(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Alternative laminate evaluation methods
function hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Array{Float64,1},xc::Array{Float64,1},yt::Array{Float64,1},
  yc::Array{Float64,1},s::Array{Float64,1},imat::Array{Float64,1})
  nlam = length(imat)
  matfail = zeros(Float64,4,nlam)
  sf = zeros(Float64,4,nlam)
  for i = 1:nlam
    matfail[:,i],sf[:,i] = hashinrotem(sigma1[i],sigma2[i],tau12[i],xt[imat[i]],
      xc[imat[i]],yt[imat[i]],yc[imat[i]],s[imat[i]])
  end
  return matfail,sf
end
hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Array{Float64,1},xc::Array{Float64,1},yt::Array{Float64,1},
  yc::Array{Float64,1},s::Array{Float64,1},lam::laminate) = hashinrotem(sigma1,
  sigma2,tau12,xt,xc,yt,yc,s,lam.imat)
hashinrotem(stress::Array{Float64,2},xt::Array{Float64,1},xc::Array{Float64,1},
  yt::Array{Float64,1},yc::Array{Float64,1},s::Array{Float64,1},
  imat::Array{Float64,1}) = hashinrotem(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,imat)
hashinrotem(stress::Array{Float64,2},xt::Array{Float64,1},xc::Array{Float64,1},
  yt::Array{Float64,1},yc::Array{Float64,1},s::Array{Float64,1},
  lam::laminate) = hashinrotem(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,lam.imat)
# Alternative laminate evaluation methods using mat::Array{material,1}
function hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},mat::material,imat::Array{Float64,1})
  nlam = length(imat)
  matfail = zeros(Float64,4,nlam)
  for i = 1:nlam
    matfail[:,i] = hashinrotem(sigma1[i],sigma2[i],tau12[i],mat[imat[i]].xt,
      mat[imat[i]].xc,mat[imat[i]].yt,mat[imat[i]].yc,mat[imat[i]].s)
  end
  return matfail
end
hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},tau12::Array{Float64,1},
  mat::material,lam::laminate) = hashinrotem(sigma1,sigma2,tau12,mat,lam.imat)
hashinrotem(stress::Array{Float64,2},mat::material,imat::Array{Float64,1}) =
  hashinrotem(stress[1,:],stress[2,:],stress[3,:],mat,imat)
hashinrotem(stress::Array{Float64,2},mat::material,lam::laminate) =
  hashinrotem(stress[1,:],stress[2,:],stress[3,:],mat,lam.imat)
