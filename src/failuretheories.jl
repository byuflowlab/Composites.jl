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
function maxstress(sigma1,sigma2,tau12,sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
  return = [sigma1[i]/sigma1tu,
           -sigma1[i]/sigma1cu,
            sigma2[i]/sigma2tu,
           -sigma2[i]/sigma2cu,
            tau12[i]/tau12u,
           -tau12[i]/tau12u];
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
function tsaiwu(sigma1,sigma2,tau12,sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
  return (sigma1[i]^2.0/(sigma1tu*sigma1cu))+
         (sigma2[i].^2.0/(sigma2tu*sigma2cu))-
         sqrt(1.0/(sigma1tu*sigma1cu)*1.0/(sigma2tu*sigma2cu))*sigma1[i]*sigma2[i]+
         (1.0/sigma1tu-1.0/sigma1cu)*sigma1[i]+
         (1.0/sigma2tu-1.0/sigma2cu)*sigma2[i]+
         tau12[i]^2.0/tau12u^2.0
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
function hashinrotem(sigma1,sigma2,tau12,sigma1tu,sigma1cu,sigma2tu,sigma2cu,tau12u)
  return [sigma1/sigma1tu,
         -sigma1/sigma1cu,
          sign(sigma2)*(sigma2^2.0/sigma2tu^2.0+tau12^2.0/tau12u^2.0),
          sign(sigma2)*-(sigma2^2.0/sigma2cu^2.0+tau12^2.0/tau12u^2.0)]
end
