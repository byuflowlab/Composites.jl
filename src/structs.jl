"""
    `Material{R<:Real}(e1,e2,g12,nu12,rho,xt,xc,yt,yc,s,t)`
Returns composite type containing composite (in-plane) material properties.
# Arguments/Fields
- `e1`: E1
- `e2`: E2
- `g12`: G12
- `nu12`: poisson's Ratio
- `rho`: density
- `xt`: longitudinal tensile ultimate strength
- `xc`: longitudinal compressive ultimate strength
- `yt`: transverse tensile ultimate strength
- `yc`: transverse compressive ultimate strength
- `s`: shear ultimate strength
- `t`: ply thickness
"""
struct Material{R<:Real}
    e1::R
    e2::R
    g12::R
    nu12::R
    rho::R
    xt::R
    xc::R
    yt::R
    yc::R
    s::R
    t::R
end

Base.convert(::Type{Material{R}}, mat::Material) where R<:Real = Material{R}(mat)
(::Type{Material{R}})(mat::Material) where R<:Real = Material{R}(
    [getfield(mat,field) for field in fieldnames(Material)]...)

"""
    `Laminate{R<:Real}(matid, nply, tply, theta)`
Returns struct containing Laminate properties
# Arguments/Fields
- `matid::Array{Int, 1}`: Material id for each lamina
- `nply::Array{Int, 1}`: number of plies in each lamina
- `tply::Array{R, 1}`: ply thickness for each lamina
- `theta::Array{Int, 1}`: orientation of each lamina (degrees)
"""
struct Laminate{R<:Real}
    matid::Array{Int,1}
    nply::Array{Int,1}
    tply::Array{R,1}
    theta::Array{R,1}
end

Base.convert(::Type{Laminate{R}}, lam::Laminate) where R<:Real = Laminate{R}(lam)
(::Type{Laminate{R}})(lam::Laminate) where R<:Real = Laminate{R}(
    [getfield(lam,field) for field in fieldnames(Laminate)]...)
