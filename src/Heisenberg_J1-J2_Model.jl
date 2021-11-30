"""
This file controls the J₁-J₂ model functions.
"""

using Revise
include("HeisenbergSpins.jl")
include("HeisenbergModel.jl")

struct J1_J2_ModelParameters
    J1_params::ModelParameters
    J2_ex::Float64
end

x_stripe_stagg( spin, site::Site2D ) = spin * cos( π * (site.xind - 1) )
y_stripe_stagg( spin, site::Site2D ) = spin * cos( π * (site.yind - 1) )

"""
Initialize a J1-J2 lattice with random spins, but fix the boundary 
in the stripe configuration
"""
function initialize_spins!( lattice_spins, latt_params, model_params::J1_J2_ModelParameters )
    x_proj = Spin3(1.,0.,0.)
    y_proj = Spin3(0.,1.,0.)

    for ydx ∈ 1:latt_params.Ly, xdx ∈ 1:latt_params.Lx
        site = site_index( Site2D(xdx, ydx), latt_params )
        lattice_spins[site] = x_stripe_stagg(x_proj, Site2D(xind, yind) ) 
        if xind > num_boundary_x_per_side && xind <= latt_params.Lx - num_boundary_x_per_side
            # Randomize the bulk 
            vector = x_stripe_stagg(x_proj, Site2D(xind, yind))
            ϕ, z = 2 * π * (-1. + 2. * rand()), -1. + 2. * rand()
            vector += model_params.initial_randomness * Spin3( cos(ϕ) * sqrt(1 - z^2), sin(ϕ) * sqrt(1 - z^2), z )
            lattice_spins[site] = copy(unit_spin3(vector))
        end
    end
    return
end

"""
Calculate the effective field for the
Heisenberg J₁-J₂ model. 
    * Here, params.J1_params.Jex > 0 is 
      antiferromagnetic.
    * neighbors is a tuple of neighbor 
      tables with the following form:
      neighbors = ( nearest, next_nearest )
"""
function effective_J1_J2_field_per_site( site, lattice_spins, params::J1_J2_ModelParameters, neighbors, one_d )
    # First include the J₁ contribution
    eff_field = effective_field_per_site( site, lattice_spins, params.J1_params, neighbors[1], one_d )
    
    # Now proceed with the J₂ contribution
    num_neighbors = length(nearest_neighbors[site, :])
    if one_d
        num_neighbors = 2
    end
    for nnn ∈ 1:num_neighbors
        eff_field += -params.J2_ex * lattice_spins[ neighbors[2][site, nnn] ]
    end
    return eff_field 
end

"""
Calculate the J1-J2 MFT spin at the site 
"""
function mft_spin_per_site( site, lattice_spins, params::J1_J2_ModelParameters, neighbors, one_d )
    eff_field = effective_J1_J2_field_per_site(site, lattice_spins, params, neighbors, one_d)
    output = avg_spin( eff_field, params.J1_params.β )
    return output
end

"""
Calculate MFT for the lattice from an initial guess
"""
function mft_lattice( lattice_spins, model_params::J1_J2_ModelParameters, latt_params::LatticeParameters, neighbors; iteration_scheme = nothing )
    new_spins = copy(lattice_spins)
    site_list = iteration_scheme
    if site_list === nothing
        site_list = 1:total_sites(latt_params)
    end
    for site ∈ site_list
        if boundary_neighbor_value != nearest_neighbors[site, 1]
            new_spins[ site ] = mft_spin_per_site( site, lattice_spins, model_params, neighbors, latt_params.Ly == 1 )  # TODO: For the J2 term, the 1d condition here is dubious.
        end
    end
    return new_spins
end