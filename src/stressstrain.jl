"""
    rotstress(stress,theta)
Rotates stress theta degrees. Stress = [sigma1,sigma2,tau12]
"""
function rotstress(stress::Array{Float64,1},theta::Float64)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0 st^2.0 2.0*st*ct; st^2.0 ct^2.0 -2.0*st*ct; -st*ct st*ct ct^2.0-st^2.0]
  return T*stress
end

"""
    rotstrain(strain,theta)
Rotates strain theta degrees. Strain = [eps1,eps2,gamma12]
"""
function rotstrain(strain::Array{Float64,1},theta::Float64)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0 st^2.0 st*ct; st^2.0 ct^2.0 -st*ct; -2.0*st*ct 2.0*st*ct ct^2.0-st^2.0]
  return T*strain
end

"""
    plystrain(tply::Array{Float64,1},nply::Array{Int64,1},theta::Array{Float64,1},
    resultantstresses::Array{Float64,1},A::Array{Float64,2},B::Array{Float64,2},D::Array{Float64,2})
Calculates material strains in each ply aligned with principal material direction
"""
function getplystrain(tply::Array{Float64,1},nply::Array{Int64,1},theta::Array{Float64,1},
  resultantstress::Array{Float64,1},A::Array{Float64,2},B::Array{Float64,2},D::Array{Float64,2})

  # number of laminas
  nlam = length(nply)

  # coordinates of top and bottom of plies
  z = getz(tply,nply)

  # stiffness matrix
  S = vcat(hcat(A,B),hcat(B',D))

  # resultant strains from resultant stresses
  resultantstrain = S\resultantstress

  # local strains at top and bottom of each lamina
  localstrain = zeros(Float64,3,nlam+1)
  for i = 1:3
    for k = 1:nlam+1
      localstrain[i,k] = resultantstrain[i]+z[k]*resultantstrain[3+i]
    end
  end

  # rotate to align with material principal axes
  lowerplystrain = zeros(Float64,3,nlam)
  upperplystrain = zeros(Float64,3,nlam)
  for k = 1:nlam
    lowerplystrain[:,k] = rotstrain(localstrain[:,k],theta[k])
    upperplystrain[:,k] = rotstrain(localstrain[:,k+1],theta[k])
  end

  return lowerplystrain,upperplystrain
end

"""
    function getplystress(plystrain::Array{Float64,2},q::Array{Float64,3})
Calculates ply stresses from ply strains and material stiffness matrix
"""
function getplystress(plystrain::Array{Float64,2},q::Array{Float64,3})
  plystress = zeros(Float64,size(plystrain)...)
  for i = 1:size(q,3)
    plystress = q[:,:,i]*plystrain[:,i]
  end
  return plystress
end
