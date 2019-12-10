var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Composites-1",
    "page": "Home",
    "title": "Composites",
    "category": "section",
    "text": "Summary: Classical Laminate Theory and Composites Material Failure CriterionAuthors: Taylor McDonnell"
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "import Pkg\nPkg.add(\"https://github.com/byuflowlab/Composites.jl\")"
},

{
    "location": "library/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library/#Public-API-1",
    "page": "Library",
    "title": "Public API",
    "category": "section",
    "text": "Modules = [Composites]\nPrivate = false\nOrder = [:function, :type]"
},

{
    "location": "library/#Composites.compliance-Tuple{AbstractArray{#s29,2} where #s29<:Real,AbstractArray{#s28,2} where #s28<:Real,AbstractArray{#s27,2} where #s27<:Real}",
    "page": "Library",
    "title": "Composites.compliance",
    "category": "method",
    "text": "compliance(A::Array{<:Real,2}, B::Array{<:Real,2}, D::Array{<:Real,2})\n\nReturns alpha, beta, and delta\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.cylinder_bending_local_buckling-NTuple{4,Any}",
    "page": "Library",
    "title": "Composites.cylinder_bending_local_buckling",
    "category": "method",
    "text": "cylinder_bending_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, R::Real, h::Real)\n\nCalculates local buckling load. R is the cylinder radius. h is the total laminate thickness.\n\nFollows approach suggested in SP-8001, which uses the theoretical axial compression cylinder buckling equation with a correction factor specific for bending cases, except assumes an infinitely long cylinder and that wave numbers are real-valued rather than integers (in order to elimate them as variables)\n\nsee An Approximate Solution and Master Curves for Buckling of Symmetrically Laminated Composite Cylinders by Michael P. Nemeth for details on the equation for the buckling of an infinitely long orthotropic cylinder\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.cylinder_uniaxial_local_buckling",
    "page": "Library",
    "title": "Composites.cylinder_uniaxial_local_buckling",
    "category": "function",
    "text": "cylinder_uniaxial_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, R::Real, h::Real, correlation::Bool=true)\n\nCalculates local buckling load. R is the cylinder radius. h is the total laminate thickness.\n\nassumes: infinitely long cylinder, simply supported, uniaxial compression, wave numbers are real-valued rather than integers (this allows them to be eliminated as variables)\n\nsee An Approximate Solution and Master Curves for Buckling of Symmetrically Laminated Composite Cylinders by Michael P. Nemeth\n\nalso see Simple Formulas and Results for Buckling Resistance and Stiffness Design of Compression-Loaded Laminated-Composite Cylinders by Michael P. Nemeth\n\nThe correlation factor from SP-8001 is also applied, but can be disabled by setting correlation=false.\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.getABD!-Tuple{AbstractArray{#s24,2} where #s24<:Real,AbstractArray{#s25,2} where #s25<:Real,AbstractArray{#s26,2} where #s26<:Real,AbstractArray{#s27,1} where #s27<:Integer,AbstractArray{#s28,1} where #s28<:Integer,AbstractArray{#s29,1} where #s29<:Real,AbstractArray{#s30,1} where #s30<:Real,AbstractArray{#s31,1} where #s31<:(AbstractArray{#s32,2} where #s32<:Real)}",
    "page": "Library",
    "title": "Composites.getABD!",
    "category": "method",
    "text": "`getABD(matid::AbstractArray{<:Integer,1},\n    nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},\n    theta::AbstractArray{<:Real,1},\n    q::AbstractArray{<:AbstractArray{<:Real,2},1})`\n\n`getABD(lam::Laminate, q::AbstractArray{<:AbstractArray{<:Real,2},1})`\n\nReturns A, B, and D matrices\n\nArguments:\n\nmatid::AbstractArray{<:Integer,1}: material id of each lamina\nnply::AbstractArray{<:Integer,1}: number of plies in each lamina\ntply::AbstractArray{<:Real,1}: thickness of a ply (m) in each lamina\ntheta::AbstractArray{<:Real,1}: orientation (deg) of each lamina\nq::AbstractArray{<:AbstractArray{<:Real,2}}: Stiffness matrix of each lamina\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.getQ",
    "page": "Library",
    "title": "Composites.getQ",
    "category": "function",
    "text": "`getQ(e1::Real, e2::Real, g12::Real, nu12::Real, theta::Real=0.0)`\n\n`getQ(mat::Material, theta::Real)`\n\n`getQ(mat::Material, lam::Laminate)`\n\nReturns Q matrix rotated theta degrees\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.getmatfail-Tuple{Real,Real,Real,Real,Real,Real,Real,Real,AbstractString}",
    "page": "Library",
    "title": "Composites.getmatfail",
    "category": "method",
    "text": "`getmatfail(sigma1::Real, sigma2::Real, tau12::Real, xt::Real,\n    xc::Real, yt::Real, yc::Real, s::Real, method::String)`\n\n`getmatfail(sigma1::Real, sigma2::Real, tau12::Real, mat::Material,\n    method::String)`\n\n`getmatfail(stress::AbstractArray{<:Real,1}, mat::Material, method::String)`\n\n`getmatfail(stress::AbstractArray{<:Real,1},xt::Real,xc::Real,yt::Real,\n    yc::Real,s::Real, method::String)`\n\n`getmatfail(stress::AbstractArray{<:AbstractArray{<:Real,1},1},\n    mat::AbstractArray{Material}, lam::Laminate, method::String)`\n\nDetermines ply failure according to specified theory. Values greater than one signify ply failure. Choose between \"maxstress\", \"tsaiwu\", and \"hashinrotem\". Safety factors are stored in the second return argument.\n\nTheory Specific Output\n\nMaximum Ply Stress:\n\nLongitudinal tension\nLongitudinal compression\nTransverse tension\nTransverse compression\nPositive shear\nNegative shear\n\nTsai-Wu\n\nFailure Criterion\n\nHashin-Rotem\n\nFiber failure in tension\nFiber failure in compression\nMatrix failure in tension\nMatrix failure in compression\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.getplystrain-Tuple{AbstractArray{#s28,1} where #s28<:Integer,AbstractArray{#s27,1} where #s27<:Real,AbstractArray{#s26,1} where #s26<:Real,AbstractArray{#s25,1} where #s25<:Real}",
    "page": "Library",
    "title": "Composites.getplystrain",
    "category": "method",
    "text": "`getplystrain(nply::AbstractArray{<:Integer,1}, tply::AbstractArray{<:Real,1},\n    theta::AbstractArray{<:Real,1}, resultantstrain::AbstractArray{<:Real,1})`\n\n`getplystrain(lam::Laminate, resultantstrain::AbstractArray{<:Real,1})`\n\nCalculates strains in each ply aligned with principal material direction.\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.getplystress-Tuple{AbstractArray{#s27,1} where #s27<:(AbstractArray{#s26,1} where #s26<:Real),AbstractArray{#s25,1} where #s25<:(AbstractArray{#s24,2} where #s24<:Real),AbstractArray{#s23,1} where #s23<:Integer}",
    "page": "Library",
    "title": "Composites.getplystress",
    "category": "method",
    "text": "`getplystress(plystrain::AbstractArray{<:AbstractArray{<:Real,1},1},\n    q::AbstractArray{<:AbstractArray{<:Real,2},1},\n    matid::AbstractArray{<:Integer,1})`\n\n`getplystress(plystrain::AbstractArray{<:AbstractArray{<:Real,1},1},\n    q::AbstractArray{<:AbstractArray{<:Real,2},1}, lam::Laminate)`\n\nCalculates ply stresses from ply strains and material stiffness matrices (q)\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.getz-Tuple{AbstractArray{#s15,1} where #s15<:Real,AbstractArray{#s16,1} where #s16<:Integer}",
    "page": "Library",
    "title": "Composites.getz",
    "category": "method",
    "text": "`getz(tply::AbstractArray{<:Real,1}, nply::AbstractArray{<:Integer,1})`\n\n`getz(lam::Laminate)`\n\nReturns a laminate\'s z-coordinates (coordinates of top and bottom of laminas) given the thickness of plies in each lamina and the number of plies in each lamina.\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.plate_shear_local_buckling-Tuple{AbstractArray{#s30,2} where #s30<:Real,AbstractArray{#s29,2} where #s29<:Real,Real}",
    "page": "Library",
    "title": "Composites.plate_shear_local_buckling",
    "category": "method",
    "text": "plate_shear_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, b::Real)\n\nCalculates local buckling load and strains. b is the panel width.\n\nassumes: balanced, symmetric, large aspect ratio, simply supported, uniaxial compression, flat rectangular plate\n\nsee chapter on Structural Component Design Techniques from Alastair Johnson section 6.2: Design of composite panels\n\nprobably assumes wave numbers are continuous in order to eliminate them\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.plate_uniaxial_local_buckling-Tuple{AbstractArray{#s30,2} where #s30<:Real,AbstractArray{#s29,2} where #s29<:Real,Real}",
    "page": "Library",
    "title": "Composites.plate_uniaxial_local_buckling",
    "category": "method",
    "text": "plate_uniaxial_local_buckling(A::Array{<:Real,2}, D::Array{<:Real,2}, b::Real)\n\nCalculates local buckling load and strains. b is the panel width.\n\nassumes: balanced, symmetric, large aspect ratio, simply supported, uniaxial compression, flat rectangular plate. buckling strain is output as a positive value\n\nsee chapter on Structural Component Design Techniques from Alastair Johnson section 6.2: Design of composite panels\n\nprobably assumes wave numbers are continuous in order to eliminate them\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.rotQ!-Tuple{AbstractArray{#s15,2} where #s15<:Real,Real}",
    "page": "Library",
    "title": "Composites.rotQ!",
    "category": "method",
    "text": "`rotQ(q::AbstractArray{<:Real,2}, theta::Real)`\n\n`rotQ(qxx::Real, q12::Real, q22::Real, q66::Real, theta::Real)`\n\nRotates Q matrix by theta degrees\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.rotstrain!-Tuple{AbstractArray{#s31,1} where #s31<:Real,Real}",
    "page": "Library",
    "title": "Composites.rotstrain!",
    "category": "method",
    "text": "`rotstrain(strain,theta)`\n\nRotates strain theta degrees. strain = [eps1,eps2,gamma12]\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.rotstress!-Tuple{AbstractArray{#s31,1} where #s31<:Real,Real}",
    "page": "Library",
    "title": "Composites.rotstress!",
    "category": "method",
    "text": "`rotstress(stress,theta)`\n\nRotates stress theta degrees. stress = [sigma1,sigma2,tau12]\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.Laminate",
    "page": "Library",
    "title": "Composites.Laminate",
    "category": "type",
    "text": "`Laminate{R<:Real}(matid, nply, tply, theta)`\n\nReturns struct containing Laminate properties\n\nArguments/Fields\n\nmatid::Array{Int, 1}: Material id for each lamina\nnply::Array{Int, 1}: number of plies in each lamina\ntply::Array{R, 1}: ply thickness for each lamina\ntheta::Array{Int, 1}: orientation of each lamina (degrees)\n\n\n\n\n\n"
},

{
    "location": "library/#Composites.Material",
    "page": "Library",
    "title": "Composites.Material",
    "category": "type",
    "text": "`Material{R<:Real}(e1,e2,g12,nu12,rho,xt,xc,yt,yc,s,t)`\n\nReturns composite type containing composite (in-plane) material properties.\n\nArguments/Fields\n\ne1: E1\ne2: E2\ng12: G12\nnu12: poisson\'s Ratio\nrho: density\nxt: longitudinal tensile ultimate strength\nxc: longitudinal compressive ultimate strength\nyt: transverse tensile ultimate strength\nyc: transverse compressive ultimate strength\ns: shear ultimate strength\nt: ply thickness\n\n\n\n\n\n"
},

{
    "location": "library/#Internal-API-1",
    "page": "Library",
    "title": "Internal API",
    "category": "section",
    "text": "Modules = [Composites]\nPublic = false\nOrder = [:function, :type]"
},

]}
