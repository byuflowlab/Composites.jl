"""
    maxstress(sigma1,sigma2,tau12,xt,xc,yt,yc,s)
    maxstress(stress,xt,xc,yt,yc,s)
    maxstress(sigma1,sigma2,tau12,mat)
    maxstress(stress,mat)
Determines ply failure according to maximum stress failure theory.

# Arguments
- `sigma1`: Ply stress in longitudinal direction (+: tension, -:compression)
- `sigma2`: Ply stress in transverse direction (+: tension, -:compression)
- `tau12`: Ply in plain shear stress
- `xt`: Ply longitudinal ultimate strength in tension
- `xc`: Ply longitudinal ultimate strength in compression
- `yt`: Ply transverse ultimate strength in tension
- `yc`: Ply transverse ultimate strength in compression
- `s`: Ply ultimate in plane shear strength
- `stress`: First three rows correspond to sigma1,sigma2,tau12 respectively
- `mat`: materials struct

Returns array of results with values greater than one corresponding to failure
arranged as follows:
1. Longitudinal tension
2. Longitudinal compression
3. Transverse tension
4. Transverse compression
5. Positive shear
6. Negative shear
"""
function maxstress(sigma1::Float64,sigma2::Float64,tau12::Float64,xt::Float64,
  xc::Float64,yt::Float64,yc::Float64,s::Float64)
  return [sigma1/xt,
         -sigma1/xc,
          sigma2/yt,
         -sigma2/yc,
          tau12/s,
         -tau12/s]
end
maxstress(stress::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = maxstress(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
maxstress(sigma1::Float64,sigma2::Float64,tau12::Float64,mat::material) =
  maxstress(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
maxstress(stress::Array{Float64,1},mat::material) = maxstress(stress[1],
  stress[2],stress[3],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

function maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64)
  return hcat(sigma1/xt,
         -sigma1/xc,
          sigma2/yt,
         -sigma2/yc,
          tau12/s,
         -tau12/s)
end
maxstress(stress::Array{Float64,2},xt::Array{Float64,1},xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = maxstress(stress[1,:],stress[2,:],stress[3,:],xt,xc,yt,yc,s)
maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},tau12::Array{Float64,1},
  mat::material) = maxstress(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
maxstress(stress::Array{Float64,2},mat::material) = maxstress(stress[1,:],
  stress[2,:],stress[3,:],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

"""
    tsaiwu(sigma1,sigma2,tau12,xt,xc,yt,yc,s)
    tsaiwu(stress,xt,xc,yt,yc,s)
    tsaiwu(sigma1,sigma2,tau12,mat)
    tsaiwu(stress,mat)
Determines ply failure according to tsai-wu failure theory. Values above one
correspond to ply failure.
# Arguments
- `sigma1`: Ply stress in longitudinal direction (+: tension, -:compression)
- `sigma2`: Ply stress in transverse direction (+: tension, -:compression)
- `tau12`: Ply in plain shear stress
- `xt`: Ply longitudinal ultimate strength in tension
- `xc`: Ply longitudinal ultimate strength in compression
- `yt`: Ply transverse ultimate strength in tension
- `yc`: Ply transverse ultimate strength in compression
- `s`: Ply ultimate in plane shear strength
- `stress`: First three rows correspond to sigma1,sigma2,tau12 respectively
- `mat`: materials struct
"""
function tsaiwu(sigma1::Float64,sigma2::Float64,tau12::Float64,xt::Float64,
  xc::Float64,yt::Float64,yc::Float64,s::Float64)
  return (sigma1^2.0/(xt*xc))+
         (sigma2^2.0/(yt*yc))-
         sqrt(1.0/(xt*xc)*1.0/(yt*yc))*sigma1*sigma2+
         (1.0/xt-1.0/xc)*sigma1+
         (1.0/yt-1.0/yc)*sigma2+
         tau12^2.0/s^2.0
end
tsaiwu(stress::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = tsaiwu(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
tsaiwu(sigma1::Float64,sigma2::Float64,tau12::Float64,mat::material) = tsaiwu(
  sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
tsaiwu(stress::Array{Float64,1},mat::material) = tsaiwu(stress[1],stress[2],
  stress[3],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

function tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64)
  return (sigma1.^2.0./(xt*xc))+
         (sigma2.^2.0./(yt*yc))-
         sqrt(1.0/(xt*xc)*1.0/(yt*yc)).*sigma1.*sigma2+
         (1.0/xt-1.0/xc).*sigma1+
         (1.0/yt-1.0/yc).*sigma2+
         tau12.^2.0/s^2.0
end
tsaiwu(stress::Array{Float64,2},xt::Array{Float64,1},xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = tsaiwu(stress[1,:],stress[2,:],stress[3,:],xt,xc,yt,yc,s)
tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},tau12::Array{Float64,1},
  mat::material) = tsaiwu(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
tsaiwu(stress::Array{Float64,2},mat::material) = tsaiwu(stress[1,:],stress[2,:],
  stress[3,:],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

"""
    hasinrotem(sigma1,sigma2,tau12,xt,xc,yt,yc,s)
    hasinrotem(stress,xt,xc,yt,yc,s)
    hasinrotem(sigma1,sigma2,tau12,mat)
    hasinrotem(stress,mat)
Determines ply failure according to Hashin-Rotem failure theory
# Arguments
- `sigma1`: Ply stress in longitudinal direction (+: tension, -:compression)
- `sigma2`: Ply stress in transverse direction (+: tension, -:compression)
- `tau12`: Ply in plain shear stress
- `xt`: Ply longitudinal ultimate strength in tension
- `xc`: Ply longitudinal ultimate strength in compression
- `yt`: Ply transverse ultimate strength in tension
- `yc`: Ply transverse ultimate strength in compression
- `s`: Ply ultimate in plane shear strength
- `stress`: First three rows correspond to sigma1,sigma2,tau12 respectively
- `mat`: materials struct

Returns array of results with values greater than one corresponding to failure
arranged as follows:
1. Fiber failure in tension
2. Fiber failure in compression
3. Matrix failure in tension
4. Matrix failure in compression
"""
function hashinrotem(sigma1::Float64,sigma2::Float64,tau12::Float64,xt::Float64,
  xc::Float64,yt::Float64,yc::Float64,s::Float64)
  return [sigma1/xt,
         -sigma1/xc,
          sign(sigma2)*(sigma2^2.0/yt^2.0+tau12^2.0/s^2.0),
          sign(sigma2)*-(sigma2^2.0/yc^2.0+tau12^2.0/s^2.0)]
end
hashinrotem(stress::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = hashinrotem(stress[1],stress[2],stress[3],xt,xc,yt,yc,s)
hashinrotem(sigma1::Float64,sigma2::Float64,tau12::Float64,mat::material) =
  hashinrotem(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
hashinrotem(stress::Array{Float64,1},mat::material) = hashinrotem(stress[1],
  stress[2],stress[3],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)

function hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},xt::Float64,xc::Float64,yt::Float64,
  yc::Float64,s::Float64)
  return hcat(sigma1./xt,
         -sigma1./xc,
          sign(sigma2)*(sigma2.^2.0/yt^2.0+tau12.^2.0./s^2.0),
          sign(sigma2)*-(sigma2.^2.0/yc^2.0+tau12.^2.0./s^2.0))
end
hashinrotem(stress::Array{Float64,2},xt::Array{Float64,1},xc::Float64,yt::Float64,
  yc::Float64,s::Float64) = hashinrotem(stress[1,:],stress[2,:],stress[3,:],xt,xc,yt,yc,s)
hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},tau12::Array{Float64,1},
  mat::material) = hashinrotem(sigma1,sigma2,tau12,mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
hashinrotem(stress::Array{Float64,2},mat::material) = hashinrotem(stress[1,:],
  stress[2,:],stress[3,:],mat.xt,mat.xc,mat.yt,mat.yc,mat.s)
