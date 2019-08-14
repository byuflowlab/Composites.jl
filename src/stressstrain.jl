"""
    `rotstress(stress,theta)`

Rotates `stress` `theta` degrees. `stress` = [sigma1,sigma2,tau12]
"""
function rotstress!(stress::AbstractArray{<:Real,1},theta::Real)
  c = cosd(theta)
  s = sind(theta)
  c2 = c^2
  s2 = s^2
  cs = c*s

  rot11 = c2
  rot12 = s2
  rot13 = 2*cs
  rot21 = s2
  rot22 = c2
  rot23 = -2*cs
  rot31 = -cs
  rot32 = cs
  rot33 = c2-s2

  tmp1 = stress[1]
  tmp2 = stress[2]
  tmp3 = stress[3]

  stress[1] = rot11*tmp1 + rot12*tmp2 + rot13*tmp3
  stress[2] = rot21*tmp1 + rot22*tmp2 + rot23*tmp3
  stress[3] = rot31*tmp1 + rot32*tmp2 + rot33*tmp3

  return stress
end

function rotstress(stress::AbstractArray{<:Real,1}, theta::Real)
  tmp = copy(stress)
  return rotstress!(tmp, theta)
end

"""
    `rotstrain(strain,theta)`

Rotates `strain` `theta` degrees. `strain` = [eps1,eps2,gamma12]
"""
function rotstrain!(strain::AbstractArray{<:Real,1},theta::Real)
  c = cosd(theta)
  s = sind(theta)
  c2 = c^2
  s2 = s^2
  cs = c*s

  rot11 = c2
  rot12 = s2
  rot13 = cs
  rot21 = s2
  rot22 = c2
  rot23 = -cs
  rot31 = -2*cs
  rot32 = 2*cs
  rot33 = c2-s2

  tmp1 = strain[1]
  tmp2 = strain[2]
  tmp3 = strain[3]

  strain[1] = rot11*tmp1 + rot12*tmp2 + rot13*tmp3
  strain[2] = rot21*tmp1 + rot22*tmp2 + rot23*tmp3
  strain[3] = rot31*tmp1 + rot32*tmp2 + rot33*tmp3

  return strain
end

function rotstrain(strain::AbstractArray{<:Real,1},theta::Real)
  tmp = copy(strain)
  return rotstrain!(tmp, theta)
end

"""
    `getplystrain(nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
        theta::AbstractArray{<:Real,1}, resultantstrain::AbstractArray{<:Real,1})`

    `getplystrain(lam::Laminate, resultantstrain::AbstractArray{<:Real,1})`

Calculates strains in each ply aligned with principal material direction.
"""
function getplystrain(nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},
    theta::AbstractArray{<:Real,1}, resultantstrain::AbstractArray{<:Real,1})

    # number of laminas
    nlam = length(nply)

    # coordinates of top and bottom of plies
    z = getz(tply, nply)

    # local strains at top and bottom of each lamina
    localstrain = [resultantstrain[1:3]+z[k]*resultantstrain[4:6] for k=1:nlam+1]

    # rotate to align with material principal axes
    lowerplystrain = [rotstrain(localstrain[k], theta[k]) for k = 1:nlam]
    upperplystrain = [rotstrain(localstrain[k+1],theta[k]) for k = 1:nlam]

  return lowerplystrain, upperplystrain
end
getplystrain(lam::Laminate, resultantstrain::AbstractArray{<:Real,1}) =
    getplystrain(lam.nply, lam.tply, lam.theta, resultantstrain)

"""
    `getplystress(plystrain::AbstractArray{<:AbstractArray{<:Real,1},1},
        q::AbstractArray{<:AbstractArray{<:Real,2},1},
        matid::AbstractArray{<:Integer,1})`

    `getplystress(plystrain::AbstractArray{<:AbstractArray{<:Real,1},1},
        q::AbstractArray{<:AbstractArray{<:Real,2},1}, lam::Laminate)`

Calculates ply stresses from ply strains and material stiffness matrices (`q`)
"""
function getplystress(plystrain::AbstractArray{<:AbstractArray{<:Real,1},1},
    q::AbstractArray{<:AbstractArray{<:Real,2},1},
    matid::AbstractArray{<:Integer,1})

    return [q[matid[i]]*plystrain[i] for i in 1:length(matid)]
end
getplystress(plystrain::AbstractArray{<:AbstractArray{<:Real,1},1},
    q::AbstractArray{<:AbstractArray{<:Real,2},1}, lam::Laminate) =
    getplystress(plystrain, q, lam.matid)
