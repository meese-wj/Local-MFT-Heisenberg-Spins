"""
This code houses the helper functions for the
Heisenberg J₁-J₂ model with magnetoelastic
coupling, a nematic anisotropy, and uniaxial
anisotropy.
"""

include("Heisenberg_J1-J2_Model.jl")

struct MagElastic_Stripe_Params
    J1J2_params::J1_J2_ModelParameters
    λ::Float64   # magnetoelastic coupling
    ε::Float64   # strength of nematic domain
    γ::Float64   # single-ion uniaxial anisotropy
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
Calculate the nematicity as a function of position.
Right now, keep it as a two steps.
"""
function nematicity( bond::Point, ε, Lx )
    center::Float64 = 1 + (Lx-1)/2.
    width::Float64  = (Lx - 1)/3.
    xmin, xmax = floor(center - 0.5 * width), ceil(center + 0.5 * width)
    # if abs( bond.xind - center ) < width/2 - 1
    # if bond.xind >= xmin && bond.xind <= xmax
    if bond.xind > num_boundary_x_per_side + 1 && bond.xind < Lx - num_boundary_x_per_side
        return -ε
    end
    return ε
end

"""
Calculate the nematic anisotropy component to the effective field.
    
    𝐡eff = -ε( 𝐒ᵢ₊ₓ - 𝐒ᵢ₊y ), 

so ε < 0 pulls the spins towards the vertical stripe state with 
ordering vector 𝐐 = (π, 0).
"""
# function nematic_anisotropy_field( site, lattice_spins, model_params::MagElastic_Stripe_Params, latt_params, nearest_neighbors )
#     eff_field =  lattice_spins[ nearest_neighbors[site, 1] ] + lattice_spins[ nearest_neighbors[site, 2] ]
#     eff_field -= lattice_spins[ nearest_neighbors[site, 3] ] + lattice_spins[ nearest_neighbors[site, 4] ]
#     return -nematicity( site_xindex(site, latt_params), model_params.ε, latt_params.Lx ) * eff_field
# end
function nematic_anisotropy_field( site, lattice_spins, model_params::MagElastic_Stripe_Params, latt_params, nearest_neighbors )
    eff_field = Spin3(0., 0., 0.)
    for nn ∈ 1:length(nearest_neighbors[site, :])
        mid = midpoint( point2d_convert(site_coords( site, latt_params )), point2d_convert(site_coords( nearest_neighbors[site, nn], latt_params )) )
        term = nematicity( mid, model_params.ε, latt_params.Lx ) * lattice_spins[ nearest_neighbors[site, nn] ]
        eff_field += (1. - 2 * (nn > 2)) * term  
    end
    return -1. * eff_field
end

"""
Calculate the magnetoelastic anisotropy component to the effective 
field.
    
    𝐡eff = -λ( Sᵢ₊ₓˣ ̂eˣ + Sᵢ₊yʸ ̂eʸ )

so λ > 0 pulls the spins towards the horizontal stripe state with 
ordering vector 𝐐 = (0, π) and spin projection ̂eˣ, or towards the 
vertical stripe state with 𝐐 = (π, 0) and spin projection ̂eʸ.
"""
function magnetoelastic_anisotropy_field( site, lattice_spins, model_params::MagElastic_Stripe_Params, nearest_neighbors )
    eff_field =  projx( lattice_spins[ nearest_neighbors[site, 1] ] ) + projx( lattice_spins[ nearest_neighbors[site, 2] ] )
    eff_field += projy( lattice_spins[ nearest_neighbors[site, 3] ] ) + projy( lattice_spins[ nearest_neighbors[site, 4] ] )
    return -model_params.λ * eff_field
end

"""
Calculate the total effective field at each site from all
contributions.
"""
function effective_4_State_Stripe_field_per_site( site, lattice_spins, params::MagElastic_Stripe_Params, latt_params, neighbors, one_d )
    # First grab the J₁-J₂ part
    eff_field = effective_J1_J2_field_per_site( site, lattice_spins, params.J1J2_params, neighbors, one_d )
    # Next grab the nematic part
    eff_field += nematic_anisotropy_field( site, lattice_spins, params, latt_params, neighbors[1] )
    # Then grab the magnetoelastic part
    eff_field += magnetoelastic_anisotropy_field( site, lattice_spins, params, neighbors[1] )
    # Finally include the on-site uniaxial anisotropy 
    eff_field += uniaxial_anisotropy_field( lattice_spins[site], params )
    return eff_field
end

"""
Calculate the J1-J2 MFT spin at the site 
"""
function mft_spin_per_site( site, lattice_spins, params::MagElastic_Stripe_Params, latt_params, neighbors, one_d )
    eff_field = effective_4_State_Stripe_field_per_site(site, lattice_spins, params, latt_params, neighbors, one_d)
    output = avg_spin( eff_field, params.J1J2_params.J1_params.β )
    return output
end

"""
Sweep the lattice and compute the total energy
"""
function mft_energy_of_system( lattice_spins, params::MagElastic_Stripe_Params, latt_params, neighbors, one_d )
    energy = 0.
    for site ∈ 1:length(lattice_spins)
        if boundary_neighbor_value != neighbors[1][site, 1]
            eff_field = effective_4_State_Stripe_field_per_site( site, lattice_spins, params, latt_params, neighbors, one_d )
            energy += mft_energy_per_spin( eff_field, lattice_spins[site] )
        end
    end
    return energy / total_bulk_sites(latt_params)
end

"""
Calculate MFT for the lattice from an initial guess
"""
function mft_lattice( lattice_spins, model_params::MagElastic_Stripe_Params, latt_params::LatticeParameters, neighbors; iteration_scheme = nothing )
    new_spins = copy(lattice_spins)
    site_list = iteration_scheme
    if site_list === nothing
        site_list = 1:total_sites(latt_params)
    end
    for site ∈ site_list
        if boundary_neighbor_value != neighbors[1][site, 1]
            new_spins[ site ] = mft_spin_per_site( site, lattice_spins, model_params, latt_params, neighbors, latt_params.Ly == 1 )  # TODO: For the J2 term, the 1d condition here is dubious.
        end
    end
    return new_spins
end
