"""
    `material{R<:Real}(e1,e2,g12,nu12,rho,xt,xc,yt,yc,s,t)`
Returns struct containing composite (in-plane) material properties.
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
struct material{R<:Real}
    e1::R #E1
    e2::R #E2
    g12::R #G12
    nu12::R #nu12
    rho::R #density
    xt::R #longitudinal tensile ultimate strength (Pa)
    xc::R #longitudinal compressive ultimate strength (Pa)
    yt::R #transverse tensile ultimate strength (Pa)
    yc::R #transverse compressive ultimate strength (Pa)
    s::R #shear ultimate strength (Pa)
    t::R # ply thickness
end

Base.convert(::Type{material{R}}, mat::material) where R<:Real = material{R}(mat)
(::Type{material{R}})(mat::material) where R<:Real = material{R}(
    [getfield(mat,field) for field in fieldnames(material)]...)

"""
    `laminate{I<:Integer,R<:Real}(matid, nply, tply, theta)`
Returns struct containing laminate properties
# Arguments/Fields
- `matid::Array{I,1}`: material id for each lamina
- `nply::Array{I,1}`: number of plies in each lamina
- `tply::Array{R,1}`: ply thickness for each lamina
- `theta::Array{I,1}`: orientation of each lamina (degrees)
"""
struct laminate{R<:Real}
    matid::Array{Int,1}
    nply::Array{Int,1}
    tply::Array{R,1}
    theta::Array{R,1}
end

Base.convert(::Type{laminate{R}}, lam::laminate) where R<:Real = laminate{R}(lam)
(::Type{laminate{R}})(lam::laminate) where R<:Real = laminate{R}(
    [getfield(lam,field) for field in fieldnames(laminate)]...)
