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
Initialize a J1-J2 lattice with random spins, but fix the boundary 
in the stripe configuration
"""
function initialize_spins!( lattice_spins, latt_params, model_params::MagElastic_Stripe_Params )
    initialize_spins!(lattice_spins, latt_params, model_params.J1J2_params)
end

"""
Calculate the uniaxial anisotropy component to the effective field.
    
    𝐡eff = -γSᶻ̂e³, so γ < 0 pulls the spins towards the z axis
"""
uniaxial_anisotropy_field( spin::Spin3, model_params::MagElastic_Stripe_Params ) = -model_params.γ * projz( spin )

"""
Calculate the nematic anisotropy component to the effective field.
    
    𝐡eff = -ε( 𝐒ᵢ₊ₓ - 𝐒ᵢ₊y ), 

so ε < 0 pulls the spins towards the horizontal stripe state with 
ordering vector 𝐐 = (0, π).
"""
function nematic_anisotropy_field( site, lattice_spins, model_params::MagElastic_Stripe_Params, nearest_neighbors )
    eff_field =  lattice_spins[ nearest_neighbors[site, 1] ] + lattice_spins[ nearest_neighbors[site, 2] ]
    eff_field -= lattice_spins[ nearest_neighbors[site, 3] ] + lattice_spins[ nearest_neighbors[site, 4] ]
    return -model_params.ε * eff_field
end

"""
Calculate the magnetoelastic anisotropy component to the effective 
field.
    
    𝐡eff = -λ( Sᵢ₊ₓˣ ̂eˣ + Sᵢ₊yʸ ̂eʸ )

so λ < 0 pulls the spins towards the horizontal stripe state with 
ordering vector 𝐐 = (0, π) and spin projection ̂eˣ, or towards the 
vertical stripe state with 𝐐 = (π, 0) and spin projection ̂eʸ.
"""
function magnetoelastic_anisotropy_field( site, lattice_spins, model_params::MagElastic_Stripe_Params, nearest_neighbors )
    eff_field =  projx( lattice_spins[ nearest_neighbors[site, 1] ] ) + projx( lattice_spins[ nearest_neighbors[site, 2] ] )
    eff_field += projy( lattice_spins[ nearest_neighbors[site, 3] ] ) + projy( lattice_spins[ nearest_neighbors[site, 4] ] )
    return -model_params.λ * eff_field
end