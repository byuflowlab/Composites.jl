"""
    rotstress(stress,theta)
Rotates stress theta degrees. Stress = [sigma1,sigma2,tau12]
"""
function rotstress(stress::Array{Real,1},theta::Real)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0 st^2.0 2.0*st*ct; st^2.0 ct^2.0 -2.0*st*ct; -st*ct st*ct ct^2.0-st^2.0]
  return T*stress
end

"""
    rotstrain(strain,theta)
Rotates strain theta degrees. Strain = [eps1,eps2,gamma12]
"""
function rotstrain(strain::Array{Real,1},theta::Real)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0 st^2.0 st*ct; st^2.0 ct^2.0 -st*ct; -2.0*st*ct 2.0*st*ct ct^2.0-st^2.0]
  return T*strain
end

"""
    plystrain(tply::Array{Real,1},nply::Array{Int64,1},theta::Array{Real,1},
    resultantstresses::Array{Real,1},A::Array{Real,2},B::Array{Real,2},D::Array{Real,2})
    getplystrain(lam::laminate,resultantstrain::Array{Real,1})
Calculates material strains in each ply aligned with principal material direction
"""
function getplystrain(nply::Array{Int64,1},tply::Array{Real,1},theta::Array{Real,1},
  resultantstrain::Array{Real,1})

  # number of laminas
  nlam = length(nply)

  # coordinates of top and bottom of plies
  z = getz(tply,nply)

  # local strains at top and bottom of each lamina
  localstrain = zeros(Real,3,nlam+1)
  for k = 1:nlam+1
    localstrain[:,k] = resultantstrain[1:3]+z[k]*resultantstrain[4:6]
  end

  # rotate to align with material principal axes
  lowerplystrain = zeros(Real,3,nlam)
  upperplystrain = zeros(Real,3,nlam)
  for k = 1:nlam
    lowerplystrain[:,k] = rotstrain(localstrain[:,k],theta[k])
    upperplystrain[:,k] = rotstrain(localstrain[:,k+1],theta[k])
  end

  return lowerplystrain,upperplystrain
end
getplystrain(lam::laminate,resultantstrain::Array{Real,1}) = getplystrain(
  lam.nply,lam.tply,lam.theta,resultantstrain)

"""
    getplystress(plystrain::Array{Real,2},q::Array{Real,3})
    getplystress(plystrain::Array{Real,2},q,lam::laminate)
Calculates ply stresses from ply strains and material stiffness matrix
"""
function getplystress(plystrain::Array{Real,2},q::Array{Real,3},matid::Array{Int64,1})
  plystress = zeros(Real,3,size(plystrain,2))
  for i = 1:length(matid)
    plystress[:,i] = q[:,:,matid[i]]*plystrain[:,i]
  end
  return plystress
end
getplystress(plystrain::Array{Real,2},q::Array{Real,3},lam::laminate) = getplystress(plystrain,q,lam.matid)
