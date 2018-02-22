# Alternative single evaluation methods
maxstress(stress::Array{Real,1},xt::Real,xc::Real,yt::Real,
  yc::Real,s::Real) = maxstress(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
maxstress(sigma1::Real,sigma2::Real,tau12::Real,mat::material) =
  maxstress(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Alternative laminate evaluation methods
function maxstress(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},xt::Array{Real,1},xc::Array{Real,1},yt::Array{Real,1},
  yc::Array{Real,1},s::Array{Real,1},imat::Array{Real,1})
  nlam = length(imat)
  matfail = zeros(Real,6,nlam)
  sf = zeros(Real,6,nlam)
  for i = 1:nlam
    matfail[:,i],sf[:,i] = maxstress(sigma1[i],sigma2[i],tau12[i],xt[imat[i]],
      xc[imat[i]],yt[imat[i]],yc[imat[i]],s[imat[i]])
  end
  return matfail,sf
end
maxstress(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},xt::Array{Real,1},xc::Array{Real,1},yt::Array{Real,1},
  yc::Array{Real,1},s::Array{Real,1},lam::laminate) = maxstress(sigma1,
  sigma2,tau12,xt,xc,yt,yc,s,lam.imat)
maxstress(stress::Array{Real,2},xt::Array{Real,1},xc::Array{Real,1},
  yt::Array{Real,1},yc::Array{Real,1},s::Array{Real,1},
  imat::Array{Real,1}) = maxstress(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,imat)
maxstress(stress::Array{Real,2},xt::Array{Real,1},xc::Array{Real,1},
  yt::Array{Real,1},yc::Array{Real,1},s::Array{Real,1},
  lam::laminate) = maxstress(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,lam.imat)
# Alternative laminate evaluation methods using mat::Array{material,1}
function maxstress(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},mat::material,imat::Array{Real,1})
  nlam = length(imat)
  matfail = zeros(Real,6,nlam)
  for i = 1:nlam
    matfail[:,i] = maxstress(sigma1[i],sigma2[i],tau12[i],mat[imat[i]].xt,
      mat[imat[i]].xc,mat[imat[i]].yt,mat[imat[i]].yc,mat[imat[i]].s)
  end
  return matfail
end
maxstress(sigma1::Array{Real,1},sigma2::Array{Real,1},tau12::Array{Real,1},
  mat::material,lam::laminate) = maxstress(sigma1,sigma2,tau12,mat,lam.imat)
maxstress(stress::Array{Real,2},mat::material,imat::Array{Real,1}) =
  maxstress(stress[1,:],stress[2,:],stress[3,:],mat,imat)
maxstress(stress::Array{Real,2},mat::material,lam::laminate) =
  maxstress(stress[1,:],stress[2,:],stress[3,:],mat,lam.imat)

# Alternative single evaluation methods
tsaiwu(stress::Array{Real,1},xt::Real,xc::Real,yt::Real,
  yc::Real,s::Real) = tsaiwu(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
tsaiwu(sigma1::Real,sigma2::Real,tau12::Real,mat::material) =
  tsaiwu(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Alternative laminate evaluation methods
function tsaiwu(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},xt::Array{Real,1},xc::Array{Real,1},yt::Array{Real,1},
  yc::Array{Real,1},s::Array{Real,1},imat::Array{Real,1})
  nlam = length(imat)
  matfail = zeros(Real,nlam)
  sf = zeros(Real,2,nlam)
  for i = 1:nlam
    matfail[i],sf[:,i] = tsaiwu(sigma1[i],sigma2[i],tau12[i],xt[imat[i]],
      xc[imat[i]],yt[imat[i]],yc[imat[i]],s[imat[i]])
  end
  return matfail,sf1,sf2
end
tsaiwu(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},xt::Array{Real,1},xc::Array{Real,1},yt::Array{Real,1},
  yc::Array{Real,1},s::Array{Real,1},lam::laminate) = tsaiwu(sigma1,
  sigma2,tau12,xt,xc,yt,yc,s,lam.imat)
tsaiwu(stress::Array{Real,2},xt::Array{Real,1},xc::Array{Real,1},
  yt::Array{Real,1},yc::Array{Real,1},s::Array{Real,1},
  imat::Array{Real,1}) = tsaiwu(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,imat)
tsaiwu(stress::Array{Real,2},xt::Array{Real,1},xc::Array{Real,1},
  yt::Array{Real,1},yc::Array{Real,1},s::Array{Real,1},
  lam::laminate) = tsaiwu(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,lam.imat)
# Alternative laminate evaluation methods using mat::Array{material,1}
function tsaiwu(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},mat::material,imat::Array{Real,1})
  nlam = length(imat)
  matfail = zeros(Real,nlam)
  for i = 1:nlam
    matfail[i] = tsaiwu(sigma1[i],sigma2[i],tau12[i],mat[imat[i]].xt,
      mat[imat[i]].xc,mat[imat[i]].yt,mat[imat[i]].yc,mat[imat[i]].s)
  end
  return matfail
end
tsaiwu(sigma1::Array{Real,1},sigma2::Array{Real,1},tau12::Array{Real,1},
  mat::material,lam::laminate) = tsaiwu(sigma1,sigma2,tau12,mat,lam.imat)
tsaiwu(stress::Array{Real,2},mat::material,imat::Array{Real,1}) =
  tsaiwu(stress[1,:],stress[2,:],stress[3,:],mat,imat)
tsaiwu(stress::Array{Real,2},mat::material,lam::laminate) =
  tsaiwu(stress[1,:],stress[2,:],stress[3,:],mat,lam.imat)


# Alternative single evaluation methods
hashinrotem(stress::Array{Real,1},xt::Real,xc::Real,yt::Real,
  yc::Real,s::Real) = hashinrotem(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
hashinrotem(sigma1::Real,sigma2::Real,tau12::Real,mat::material) =
  hashinrotem(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Alternative laminate evaluation methods
function hashinrotem(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},xt::Array{Real,1},xc::Array{Real,1},yt::Array{Real,1},
  yc::Array{Real,1},s::Array{Real,1},imat::Array{Real,1})
  nlam = length(imat)
  matfail = zeros(Real,4,nlam)
  sf = zeros(Real,4,nlam)
  for i = 1:nlam
    matfail[:,i],sf[:,i] = hashinrotem(sigma1[i],sigma2[i],tau12[i],xt[imat[i]],
      xc[imat[i]],yt[imat[i]],yc[imat[i]],s[imat[i]])
  end
  return matfail,sf
end
hashinrotem(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},xt::Array{Real,1},xc::Array{Real,1},yt::Array{Real,1},
  yc::Array{Real,1},s::Array{Real,1},lam::laminate) = hashinrotem(sigma1,
  sigma2,tau12,xt,xc,yt,yc,s,lam.imat)
hashinrotem(stress::Array{Real,2},xt::Array{Real,1},xc::Array{Real,1},
  yt::Array{Real,1},yc::Array{Real,1},s::Array{Real,1},
  imat::Array{Real,1}) = hashinrotem(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,imat)
hashinrotem(stress::Array{Real,2},xt::Array{Real,1},xc::Array{Real,1},
  yt::Array{Real,1},yc::Array{Real,1},s::Array{Real,1},
  lam::laminate) = hashinrotem(stress[1,:],stress[2,:],stress[3,:],xt,xc,
  yt,yc,s,lam.imat)
# Alternative laminate evaluation methods using mat::Array{material,1}
function hashinrotem(sigma1::Array{Real,1},sigma2::Array{Real,1},
  tau12::Array{Real,1},mat::material,imat::Array{Real,1})
  nlam = length(imat)
  matfail = zeros(Real,4,nlam)
  for i = 1:nlam
    matfail[:,i] = hashinrotem(sigma1[i],sigma2[i],tau12[i],mat[imat[i]].xt,
      mat[imat[i]].xc,mat[imat[i]].yt,mat[imat[i]].yc,mat[imat[i]].s)
  end
  return matfail
end
hashinrotem(sigma1::Array{Real,1},sigma2::Array{Real,1},tau12::Array{Real,1},
  mat::material,lam::laminate) = hashinrotem(sigma1,sigma2,tau12,mat,lam.imat)
hashinrotem(stress::Array{Real,2},mat::material,imat::Array{Real,1}) =
  hashinrotem(stress[1,:],stress[2,:],stress[3,:],mat,imat)
hashinrotem(stress::Array{Real,2},mat::material,lam::laminate) =
  hashinrotem(stress[1,:],stress[2,:],stress[3,:],mat,lam.imat)
