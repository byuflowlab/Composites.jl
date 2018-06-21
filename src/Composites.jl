"""
This module uses classical laminate theory to calculate stresses and strains in
a thin plate laminate.  Analysis is based on chapter 3 of "Design and Analysis
of Composite Structures" by Kassapoglou.
"""
module Composites

__precompile__()

# Some modules structs
include("structs.jl")
# ABD matrix calculations
include("abd.jl")
# Stress and strain calculations
include("stressstrain.jl")
# Material Failure Theories
include("failuretheories.jl")
# Buckling related functions
include("buckling.jl")

end #module
