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
    xc::Real, yt::Real, yc::Real, s::Real, method::AbstractString)

    if method == "maxstress"
        return maxstress(sigma1, sigma2, tau12, xt, xc, yt, yc, s)
    elseif method == "tsaiwu"
        return tsaiwu(sigma1, sigma2, tau12, xt, xc, yt, yc, s)
    elseif method == "hashinrotem"
        return hashinrotem(sigma1, sigma2, tau12, xt, xc, yt, yc, s)
    else
        error("Chosen method ($method) not implemented. Choose between `maxstress`,
            `tsaiwu`, and `hashinrotem`")
    end
end

getmatfail(stress::AbstractArray{<:Real,1}, mat::Material, method::AbstractString) =
    getmatfail(stress[1], stress[2], stress[3], mat.xt, mat.xc, mat.yt, mat.yc,
    mat.s, method)
getmatfail(stress::AbstractArray{<:Real,1}, xt::Real, xc::Real, yt::Real,
    yc::Real, s::Real, method::AbstractString) = getmatfail(stress[1], stress[2],
    stress[3], xt, xc, yt, yc, s, method)
getmatfail(sigma1::Real, sigma2::Real, tau12::Real, mat::Material, method::AbstractString) =
    getmatfail(sigma1, sigma2, tau12, mat.xt, mat.xc, mat.yt, mat.yc, mat.s, method)

function getmatfail(stress::AbstractArray{<:AbstractArray{<:Real,1},1},
    mat::AbstractArray{<:Material,1}, lam::Laminate, method::AbstractString)

    result = [getmatfail(stress[i][1], stress[i][2], stress[i][3],
        mat[lam.matid[i]].xt, mat[lam.matid[i]].xc,
        mat[lam.matid[i]].yt, mat[lam.matid[i]].yc,
        mat[lam.matid[i]].s, method) for i in 1:length(lam.matid)]

    matfail = [result[i][1] for i in 1:length(lam.matid)]
    sf = [result[i][2] for i in 1:length(lam.matid)]

    return matfail, sf
end

function maxstress(sigma1::Real, sigma2::Real, tau12::Real, xt::Real,
    xc::Real, yt::Real, yc::Real, s::Real)

    fail = [sigma1/xt, -sigma1/xc, sigma2/yt, -sigma2/yc, tau12/s, -tau12/s]

    safetyfactor = 1.0./fail

    return fail, safetyfactor
end
maxstress(stress::AbstractArray{<:Real,1}, mat::Material) =
    maxstress(stress[1], stress[2], stress[3], mat.xt, mat.xc, mat.yt, mat.yc, mat.s)
maxstress(stress::AbstractArray{<:Real,1}, xt::Real, xc::Real, yt::Real,
    yc::Real, s::Real) = maxstress(stress[1], stress[2], stress[3], xt, xc, yt, yc, s)
maxstress(sigma1::Real, sigma2::Real, tau12::Real, mat::Material) =
    maxstress(sigma1, sigma2, tau12, mat.xt, mat.xc, mat.yt, mat.yc, mat.s)

function maxstress(stress::AbstractArray{<:AbstractArray{<:Real,1},1},
    mat::AbstractArray{<:Material,1}, lam::Laminate)

    result = [maxstress(stress[i][1], stress[i][2], stress[i][3],
        mat[lam.matid[i]].xt, mat[lam.matid[i]].xc,
        mat[lam.matid[i]].yt, mat[lam.matid[i]].yc,
        mat[lam.matid[i]].s, method) for i in 1:length(lam.matid)]

    matfail = [result[i][1] for i in 1:length(lam.matid)]
    sf = [result[i][2] for i in 1:length(lam.matid)]

    return matfail, sf
end

function tsaiwu(sigma1::Real, sigma2::Real, tau12::Real, xt::Real,
    xc::Real, yt::Real, yc::Real, s::Real)

    a = (sigma1^2.0/(xt*xc))+(sigma2^2.0/(yt*yc))-
         sqrt(1.0/(xt*xc)*1.0/(yt*yc))*sigma1*sigma2+tau12^2.0/s^2.0
    b = (1.0/xt-1.0/xc)*sigma1+(1.0/yt-1.0/yc)*sigma2
    c = -1
    fail = a+b
    safetyfactor = abs(fail) > 1e-20 ? (-b+sqrt(b^2-4*a*c))/(2a) : Inf

    return fail, safetyfactor
end
tsaiwu(stress::AbstractArray{<:Real,1}, mat::Material) =
    tsaiwu(stress[1], stress[2], stress[3], mat.xt, mat.xc, mat.yt, mat.yc, mat.s)
tsaiwu(stress::AbstractArray{<:Real,1}, xt::Real, xc::Real, yt::Real,
    yc::Real, s::Real) = tsaiwu(stress[1], stress[2], stress[3], xt, xc, yt, yc, s)
tsaiwu(sigma1::Real, sigma2::Real, tau12::Real, mat::Material) =
    tsaiwu(sigma1, sigma2, tau12, mat.xt, mat.xc, mat.yt, mat.yc, mat.s)

function tsaiwu(stress::AbstractArray{<:AbstractArray{<:Real,1},1},
    mat::AbstractArray{<:Material,1}, lam::Laminate)

    result = [tsaiwu(stress[i][1], stress[i][2], stress[i][3],
        mat[lam.matid[i]].xt, mat[lam.matid[i]].xc,
        mat[lam.matid[i]].yt, mat[lam.matid[i]].yc,
        mat[lam.matid[i]].s, method) for i in 1:length(lam.matid)]

    matfail = [result[i][1] for i in 1:length(lam.matid)]
    sf = [result[i][2] for i in 1:length(lam.matid)]

    return matfail, sf
end

function hashinrotem(sigma1::Real, sigma2::Real, tau12::Real, xt::Real,
    xc::Real, yt::Real, yc::Real, s::Real)

    tensionfail = (sigma2^2.0/yt^2.0+tau12^2.0/s^2.0)
    compressionfail = (sigma2^2.0/yc^2.0+tau12^2.0/s^2.0)
    fail = [sigma1/xt, -sigma1/xc, sign(sigma2)*tensionfail, sign(sigma2)*-compressionfail]
    safetyfactor = [xt/sigma1, -xc/sigma1, 1/sqrt(tensionfail), sign(sigma2)*-1/sqrt(compressionfail)]

    return fail, safetyfactor
end
hashinrotem(stress::AbstractArray{<:Real,1}, mat::Material) =
    hashinrotem(stress[1], stress[2], stress[3], mat.xt, mat.xc, mat.yt, mat.yc, mat.s)
hashinrotem(stress::AbstractArray{<:Real,1}, xt::Real, xc::Real, yt::Real,
    yc::Real, s::Real) = hashinrotem(stress[1], stress[2], stress[3], xt, xc, yt, yc, s)
hashinrotem(sigma1::Real, sigma2::Real, tau12::Real, mat::Material) =
    hashinrotem(sigma1, sigma2, tau12, mat.xt, mat.xc, mat.yt, mat.yc, mat.s)

function hashinrotem(stress::AbstractArray{<:AbstractArray{<:Real,1},1},
    mat::AbstractArray{<:Material,1}, lam::Laminate)

    result = [hashinrotem(stress[i][1], stress[i][2], stress[i][3],
        mat[lam.matid[i]].xt, mat[lam.matid[i]].xc,
        mat[lam.matid[i]].yt, mat[lam.matid[i]].yc,
        mat[lam.matid[i]].s, method) for i in 1:length(lam.matid)]

    matfail = [result[i][1] for i in 1:length(lam.matid)]
    sf = [result[i][2] for i in 1:length(lam.matid)]

    return matfail, sf
end
