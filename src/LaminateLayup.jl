"""
This module uses classical laminate theory to calculate stresses and strains in
a thin plate laminate.  Analysis is based on chapter 3 of "Design and Analysis
of Composite Structures" by Kassapoglou.
"""
module LaminateLayup
export stiffness,compliance,plystrain,plystress

# ABD matrix calculations
include("abd.jl")
# Stress and strain calculations
include("stressstrain.jl")
# Material Failure Theories
include("failuretheories.jl")

end #module
