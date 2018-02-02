"""
    maxstress(stress::Array{Float64,1},mat::material)
    maxstress(stress::Array{Float64,2},mat::Array{material,1},lam::laminate)
Determines ply failure according to maximum stress failure theory.
# Arguments
- `stress`: First three rows correspond to sigma1,sigma2,tau12 respectively
- `mat`: materials struct(s)
- `lam`: laminate struct

Returns array of results with values greater than one corresponding to failure
arranged as follows:
Row 1: Longitudinal tension
Row 2: Longitudinal compression
Row 3: Transverse tension
Row 4: Transverse compression
Row 5: Positive shear
Row 6: Negative shear
"""
# Single evaluation method
function maxstress(sigma1::Float64,sigma2::Float64,tau12::Float64,xt::Float64,
  xc::Float64,yt::Float64,yc::Float64,s::Float64)
  fail = [sigma1/xt,
         -sigma1/xc,
          sigma2/yt,
         -sigma2/yc,
          tau12/s,
         -tau12/s]
  safetyfactor = 1./fail
  return fail,safetyfactor
end
maxstress(stress::Array{Float64,1},mat::material) = maxstress(stress[1],
  stress[2],stress[3],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Laminate evaluation method
function maxstress(stress::Array{Float64,2},mat::Array{material},lam::laminate)
  nlam = length(lam.matid)
  matfail = zeros(Float64,6,nlam)
  sf = zeros(Float64,6,nlam)
  for i = 1:nlam
    matfail[:,i],sf[:,i] = maxstress(stress[1,i],stress[2,i],stress[3,i],mat[lam.matid[i]].xt,
      mat[lam.matid[i]].xc,mat[lam.matid[i]].yt,mat[lam.matid[i]].yc,mat[lam.matid[i]].s)
  end
  return matfail,sf
end

"""
  tsaiwu(stress::Array{Float64,1},mat::material)
  tsaiwu(stress::Array{Float64,2},mat::Array{material,1},lam::laminate)
Determines ply failure according to tsai-wu failure theory. Values above one
correspond to ply failure.
# Arguments
- `stress`: First three rows correspond to sigma1,sigma2,tau12 respectively
- `mat`: materials struct(s)
"""
function tsaiwu(sigma1::Float64,sigma2::Float64,tau12::Float64,xt::Float64,
  xc::Float64,yt::Float64,yc::Float64,s::Float64)
  a = (sigma1^2.0/(xt*xc))+(sigma2^2.0/(yt*yc))-
    sqrt(1.0/(xt*xc)*1.0/(yt*yc))*sigma1*sigma2+tau12^2.0/s^2.0
  b = (1.0/xt-1.0/xc)*sigma1+(1.0/yt-1.0/yc)*sigma2
  c = -1
  safetyfactor = [(-b+sqrt(b^2-4*a*c))/2a,(-b-sqrt(b^2-4*a*c))/2a]
  return a+b,safetyfactor
end
tsaiwu(stress::Array{Float64,1},mat::material) = tsaiwu(stress[1],
  stress[2],stress[3],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Laminate evaluation method
function tsaiwu(stress::Array{Float64,2},mat::Array{material},lam::laminate)
  nlam = length(lam.matid)
  matfail = zeros(Float64,nlam)
  sf = zeros(Float64,2,nlam)
  for i = 1:nlam
    matfail[i],sf[:,i] = tsaiwu(stress[1,i],stress[2,i],stress[3,i],mat[lam.matid[i]].xt,
      mat[lam.matid[i]].xc,mat[lam.matid[i]].yt,mat[lam.matid[i]].yc,mat[lam.matid[i]].s)
  end
  return matfail,sf
end


"""
    hashinrotem(stress::Array{Float64,1},mat::material)
    hashinrotem(stress::Array{Float64,2},mat::Array{material,1},lam::laminate)
Determines ply failure according to Hashin-Rotem failure theory.
# Arguments
- `stress`: First three rows correspond to sigma1,sigma2,tau12 respectively
- `mat`: materials struct(s)
- `lam`: laminate struct

Returns array of results with values greater than one corresponding to failure
arranged as follows:
1. Fiber failure in tension
2. Fiber failure in compression
3. Matrix failure in tension
4. Matrix failure in compression
"""
# Single evaluation method
function hashinrotem(sigma1::Float64,sigma2::Float64,tau12::Float64,xt::Float64,
  xc::Float64,yt::Float64,yc::Float64,s::Float64)
  tensionfail = (sigma2^2.0/yt^2.0+tau12^2.0/s^2.0)
  compressionfail = (sigma2^2.0/yc^2.0+tau12^2.0/s^2.0)
  matfail = [sigma1/xt,
            -sigma1/xc,
             sign(sigma2)*tensionfail,
             sign(sigma2)*-compressionfail]
  safetyfactor = [xt/sigma1,
                  -xc/sigma1,
                  1/sqrt(tensionfail),
                  sign(sigma2)*-1/sqrt(compressionfail)]
  return matfail,safetyfactor
end
hashinrotem(stress::Array{Float64,1},mat::material) = hashinrotem(stress[1],
  stress[2],stress[3],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

# Laminate evaluation method
function hashinrotem(stress::Array{Float64,2},mat::Array{material},lam::laminate)
  nlam = length(lam.matid)
  matfail = zeros(Float64,4,nlam)
  sf = zeros(Float64,4,nlam)
  for i = 1:nlam
    matfail[:,i],sf[:,i] = hashinrotem(stress[1,i],stress[2,i],stress[3,i],mat[lam.matid[i]].xt,
      mat[lam.matid[i]].xc,mat[lam.matid[i]].yt,mat[lam.matid[i]].yc,mat[lam.matid[i]].s)
  end
  return matfail,sf
end
