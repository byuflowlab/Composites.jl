"""
    maxstress(sigma1,sigma2,tau12,sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
Determines ply failure according to maximum stress failure theory
# Arguments
- `sigma1`: Ply stress in longitudinal direction (+: tension, -:compression)
- `sigma2`: Ply stress in transverse direction (+: tension, -:compression)
- `tau12`: Ply in plain shear stress
- `sigma1tu`: Ply longitudinal ultimate strength in tension
- `sigma1cu`: Ply longitudinal ultimate strength in compression
- `sigma2tu`: Ply transverse ultimate strength in tension
- `sigma2cu`: Ply transverse ultimate strength in compression
- `tau12u`: Ply ultimate in plane shear strength

Returns array of results with values greater than one corresponding to failure
arranged as follows:
1. Longitudinal tension
2. Longitudinal compression
3. Transverse tension
4. Transverse compression
5. Positive shear
6. Negative shear
"""
function maxstress(sigma1::Float64,sigma2::Float64,tau12::Float64,sigma1tu::Float64,
  sigma1cu::Float64,sigma2tu::Float64,sigma2cu::Float64,tau12u::Float64)
  return [sigma1/sigma1tu,
         -sigma1/sigma1cu,
          sigma2/sigma2tu,
         -sigma2/sigma2cu,
          tau12/tau12u,
         -tau12/tau12u]
end

function maxstress(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},sigma1tu::Float64,sigma1cu::Float64,sigma2tu::Float64,
  sigma2cu::Float64,tau12u::Float64)
  return hcat(sigma1/sigma1tu,
         -sigma1/sigma1cu,
          sigma2/sigma2tu,
         -sigma2/sigma2cu,
          tau12/tau12u,
         -tau12/tau12u)
end

"""
    tsaiwu(sigma1,sigma2,tau12,sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
Determines ply failure according to tsai-wu failure theory. Values above one
correspond to ply failure.
# Arguments
- `sigma1`: Ply stress in longitudinal direction (+: tension, -:compression)
- `sigma2`: Ply stress in transverse direction (+: tension, -:compression)
- `tau12`: Ply in plain shear stress
- `sigma1tu`: Ply longitudinal ultimate strength in tension
- `sigma1cu`: Ply longitudinal ultimate strength in compression
- `sigma2tu`: Ply transverse ultimate strength in tension
- `sigma2cu`: Ply transverse ultimate strength in compression
- `tau12u`: Ply ultimate in plane shear strength
"""
function tsaiwu(sigma1::Float64,sigma2::Float64,tau12::Float64,sigma1tu::Float64,
  sigma1cu::Float64,sigma2tu::Float64,sigma2cu::Float64,tau12u::Float64)
  return (sigma1^2.0/(sigma1tu*sigma1cu))+
         (sigma2^2.0/(sigma2tu*sigma2cu))-
         sqrt(1.0/(sigma1tu*sigma1cu)*1.0/(sigma2tu*sigma2cu))*sigma1*sigma2+
         (1.0/sigma1tu-1.0/sigma1cu)*sigma1+
         (1.0/sigma2tu-1.0/sigma2cu)*sigma2+
         tau12^2.0/tau12u^2.0
end

function tsaiwu(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},sigma1tu::Float64,sigma1cu::Float64,sigma2tu::Float64,
  sigma2cu::Float64,tau12u::Float64)
  return (sigma1.^2.0./(sigma1tu*sigma1cu))+
         (sigma2.^2.0./(sigma2tu*sigma2cu))-
         sqrt(1.0/(sigma1tu*sigma1cu)*1.0/(sigma2tu*sigma2cu)).*sigma1.*sigma2+
         (1.0/sigma1tu-1.0/sigma1cu).*sigma1+
         (1.0/sigma2tu-1.0/sigma2cu).*sigma2+
         tau12.^2.0/tau12u^2.0
end

"""
    maxstress(sigma1,sigma2,tau12,sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
Determines ply failure according to Hashin-Rotem failure theory
# Arguments
- `sigma1`: Ply stress in longitudinal direction (+: tension, -:compression)
- `sigma2`: Ply stress in transverse direction (+: tension, -:compression)
- `tau12`: Ply in plain shear stress
- `sigma1tu`: Ply longitudinal ultimate strength in tension
- `sigma1cu`: Ply longitudinal ultimate strength in compression
- `sigma2tu`: Ply transverse ultimate strength in tension
- `sigma2cu`: Ply transverse ultimate strength in compression
- `tau12u`: Ply ultimate in plane shear strength

Returns array of results with values greater than one corresponding to failure
arranged as follows:
1. Fiber failure in tension
2. Fiber failure in compression
3. Matrix failure in tension
4. Matrix failure in compression
"""
function hashinrotem(sigma1::Float64,sigma2::Float64,tau12::Float64,sigma1tu::Float64,
  sigma1cu::Float64,sigma2tu::Float64,sigma2cu::Float64,tau12u::Float64)
  return [sigma1/sigma1tu,
         -sigma1/sigma1cu,
          sign(sigma2)*(sigma2^2.0/sigma2tu^2.0+tau12^2.0/tau12u^2.0),
          sign(sigma2)*-(sigma2^2.0/sigma2cu^2.0+tau12^2.0/tau12u^2.0)]
end

function hashinrotem(sigma1::Array{Float64,1},sigma2::Array{Float64,1},
  tau12::Array{Float64,1},sigma1tu::Float64,sigma1cu::Float64,sigma2tu::Float64,
  sigma2cu::Float64,tau12u::Float64)
  return hcat(sigma1./sigma1tu,
         -sigma1./sigma1cu,
          sign(sigma2)*(sigma2.^2.0/sigma2tu^2.0+tau12.^2.0./tau12u^2.0),
          sign(sigma2)*-(sigma2.^2.0/sigma2cu^2.0+tau12.^2.0./tau12u^2.0))
end
