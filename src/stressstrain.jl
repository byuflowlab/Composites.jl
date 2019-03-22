"""
    `rotstress(stress,theta)`

Rotates `stress` `theta` degrees. `stress` = [sigma1,sigma2,tau12]
"""
function rotstress(stress::AbstractArray{<:Real,1},theta::Real)
  ct = cosd(theta)
  st = sind(theta)
  T = [ct^2.0   st^2.0    2.0*st*ct  ;
       st^2.0   ct^2.0   -2.0*st*ct  ;
      -st*ct    st*ct   ct^2.0-st^2.0]
  return T*stress
end

"""
    `rotstrain(strain,theta)`

Rotates `strain` `theta` degrees. `strain` = [eps1,eps2,gamma12]
"""
function rotstrain(strain::AbstractArray{<:Real,1},theta::Real)
  ct = cosd(theta)
  st = sind(theta)
  T = [  ct^2.0       st^2.0         st*ct    ;
         st^2.0       ct^2.0        -st*ct    ;
       -2.0*st*ct    2.0*st*ct   ct^2.0-st^2.0]
  return T*strain
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
