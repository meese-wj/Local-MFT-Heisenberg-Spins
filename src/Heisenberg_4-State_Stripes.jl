"""
This code houses the helper functions for the
Heisenberg J₁-J₂ model with magnetoelastic
coupling, a nematic anisotropy, and uniaxial
anisotropy.
"""

include("Heisenberg_J1-J2_Model.jl")

struct MagElastic_Stripe_Params
    J1J2_params::J1_J2_ModelParameters
    λ::Float64
    ε::Float64
    γ::Float64
end

"""
Calculate
"""

