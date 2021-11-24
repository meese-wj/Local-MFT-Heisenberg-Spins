include("HeisenbergSpins.jl")
include("LangevinFunction.jl")
include("LatticeSetup.jl")

using Random

struct ModelParameters
    Jex::Float64
    β::Float64
    initial_randomness::Float64
end

temperature(params::ModelParameters) = 1/params.β

"""
Initialize lattice with random spins, but fix the boundary
"""
function initialize_spins!( lattice_spins, latt_params, model_params )
    left_boundary =  Spin3(0.,0.,1.)
    right_boundary = -1. * left_boundary
    for ydx ∈ 1:latt_params.Ly, xdx ∈ 1:latt_params.Lx
        site = site_index( Site2D(xdx, ydx), latt_params )
        if xdx <= num_boundary_x_per_side
            lattice_spins[site] = copy( left_boundary )
        elseif  xdx > latt_params.Lx - num_boundary_x_per_side
            lattice_spins[site] = copy( right_boundary )
        else
            # Randomize the bulk 
            vector = Spin3( sin( π*(xdx - 2.)/(latt_params.Lx - 2.) ), 0., cos( π*(xdx - 2.)/(latt_params.Lx - 2.) ) )
            vector += model_params.initial_randomness * Spin3( -1. + 2. * rand(), -1. + 2. * rand(), -1. + 2. * rand() )
            lattice_spins[site] = copy(unit_spin3(vector))
        end
    end
    return
end

"""
Calculate the effective field for the
Heisenberg model with nearest neighbor interactions
only. Here, Jex > 0 is antiferromagnetic.
"""
function effective_field_per_site( site, lattice_spins, params::ModelParameters, nearest_neighbors, one_d )
    eff_field = Spin3(0., 0., 0.)
    num_neighbors = length(nearest_neighbors[site, :])
    if one_d
        num_neighbors = 2
    end
    for nn ∈ 1:num_neighbors
        eff_field += -params.Jex * lattice_spins[ nearest_neighbors[site, nn] ]
    end
    return eff_field 
end

"""
Calculate the MFT spin at the site 
"""
function mft_spin_per_site( site, lattice_spins, params::ModelParameters, nearest_neighbors, one_d )
    eff_field = effective_field_per_site(site, lattice_spins, params, nearest_neighbors, one_d)
    output = unit_spin3( eff_field ) * LangevinFunction( abs(eff_field), params.β )
    return output
end
    

"""
Calculate MFT for the lattice from an initial guess
"""
function mft_lattice( lattice_spins, model_params::ModelParameters, latt_params::LatticeParameters, nearest_neighbors )
    new_spins = copy(lattice_spins)
    for site ∈ 1:total_sites(latt_params)
        if boundary_neighbor_value != nearest_neighbors[site, 1]
            new_spins[ site ] = mft_spin_per_site( site, lattice_spins, model_params, nearest_neighbors, latt_params.Ly == 1 )
        end
    end
    return new_spins
end

"""
Calculate the difference between two field configurations
"""
function average_spin_difference( field1, field2 )
    diff = field1 .- field2
    error = 0.
    for spin ∈ diff
        error += abs2(spin)
    end
    return sqrt(error) / length(diff)
end