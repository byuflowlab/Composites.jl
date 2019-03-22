"""
    `getmatfail(sigma1::Real, sigma2::Real, tau12::Real, xt::Real,
        xc::Real, yt::Real, yc::Real, s::Real, method::String)`

    `getmatfail(sigma1::Real, sigma2::Real, tau12::Real, mat::Material,
        method::String)`

    `getmatfail(stress::AbstractArray{<:Real,1}, mat::Material, method::String)`

    `getmatfail(stress::AbstractArray{<:Real,1},xt::Real,xc::Real,yt::Real,
        yc::Real,s::Real, method::String)`

    `getmatfail(stress::AbstractArray{<:AbstractArray{<:Real,1},1},
        mat::AbstractArray{Material}, lam::Laminate, method::String)`

Determines ply failure according to specified theory. Values greater than one
signify ply failure. Choose between "maxstress", "tsaiwu", and "hashinrotem".
Safety factors are stored in the second return argument.

# Theory Specific Output
Maximum Ply Stress:
1. Longitudinal tension
2. Longitudinal compression
3. Transverse tension
4. Transverse compression
5. Positive shear
6. Negative shear

Tsai-Wu
1. Failure Criterion

Hashin-Rotem
1. Fiber failure in tension
2. Fiber failure in compression
3. Matrix failure in tension
4. Matrix failure in compression
"""
function getmatfail(sigma1::Real, sigma2::Real, tau12::Real, xt::Real,
    xc::Real, yt::Real, yc::Real, s::Real, method::String)

    if method == "maxstress"
        fail = [sigma1/xt, -sigma1/xc, sigma2/yt, -sigma2/yc, tau12/s, -tau12/s]
        safetyfactor = 1.0./fail
    elseif method == "tsaiwu"
        a = (sigma1^2.0/(xt*xc))+(sigma2^2.0/(yt*yc))-
             sqrt(1.0/(xt*xc)*1.0/(yt*yc))*sigma1*sigma2+tau12^2.0/s^2.0
        b = (1.0/xt-1.0/xc)*sigma1+(1.0/yt-1.0/yc)*sigma2
        c = -1
        fail = a+b
        safetyfactor = abs(fail) > 1e-20 ? (-b+sqrt(b^2-4*a*c))/(2a) : 999999.0
    elseif method == "hashinrotem"
        tensionfail = (sigma2^2.0/yt^2.0+tau12^2.0/s^2.0)
        compressionfail = (sigma2^2.0/yc^2.0+tau12^2.0/s^2.0)
        fail = [sigma1/xt, -sigma1/xc, sign(sigma2)*tensionfail,
            sign(sigma2)*-compressionfail]
        safetyfactor = [xt/sigma1, -xc/sigma1, 1/sqrt(tensionfail),
            sign(sigma2)*-1/sqrt(compressionfail)]
    else
        error("Chosen method not implemented. Choose between `maxstress`,
            `tsaiwu`, and `hashinrotem`")
    end

    return fail, safetyfactor
end

getmatfail(stress::AbstractArray{<:Real,1}, mat::Material, method::String) =
    getmatfail(stress[1], stress[2], stress[3], mat.xt, mat.xc, mat.yt, mat.yc,
    mat.s, method)

getmatfail(stress::AbstractArray{<:Real,1}, xt::Real, xc::Real, yt::Real,
    yc::Real, s::Real, method::String) = getmatfail(stress[1], stress[2],
    stress[3], xt, xc, yt, yc, s, method)

getmatfail(sigma1::Real, sigma2::Real, tau12::Real, mat::Material, method::String) =
    getmatfail(sigma1, sigma2, tau12, mat.xt, mat.xc, mat.yt, mat.yc, mat.s, method)

function getmatfail(stress::AbstractArray{<:AbstractArray{<:Real,1},1},
    mat::AbstractArray{<:Material,1}, lam::Laminate, method::String)

    result = [getmatfail(stress[i][1], stress[i][2], stress[i][3],
        mat[lam.matid[i]].xt, mat[lam.matid[i]].xc,
        mat[lam.matid[i]].yt, mat[lam.matid[i]].yc,
        mat[lam.matid[i]].s, method) for i in 1:length(lam.matid)]

    matfail = [result[i][1] for i in 1:length(lam.matid)]
    sf = [result[i][2] for i in 1:length(lam.matid)]

    return matfail, sf
end
