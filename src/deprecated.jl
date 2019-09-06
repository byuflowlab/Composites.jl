struct material{R<:Real}
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

getQ(mat::material, theta::Real=0.0) = getQ(mat.e1, mat.e2, mat.g12, mat.nu12,
    theta)

getQ(mat::Array{<:material,1}, lam::Laminate) = getQ.(mat, lam.theta)
