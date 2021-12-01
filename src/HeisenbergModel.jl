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
checkerboard_stagger(Jex, xdx, ydx) = ( ( -1. * Jex / abs(Jex) ) ^ ( (xdx - 1) + (ydx - 1) ) )::Float64

"""
Initialize lattice with random spins, but fix the boundary
"""
function initialize_spins!( lattice_spins, latt_params, model_params )
    left_boundary =  Spin3(0.,0.,1.)
    right_boundary = -1. * left_boundary
    for ydx ∈ 1:latt_params.Ly, xdx ∈ 1:latt_params.Lx
        site = site_index( Site2D(xdx, ydx), latt_params )
        if xdx <= num_boundary_x_per_side
            lattice_spins[site] = left_boundary  * checkerboard_stagger(model_params.Jex, xdx, ydx)
        elseif  xdx > latt_params.Lx - num_boundary_x_per_side
            lattice_spins[site] = right_boundary * checkerboard_stagger(model_params.Jex, xdx, ydx)
        else
            # Randomize the bulk 
            vector = Spin3( 0., sin( π*(xdx - 2.)/(latt_params.Lx - 1.) ), cos( π*(xdx - 2.)/(latt_params.Lx - 1.) ) )
            # @show ( -1. * model_params.Jex / abs(model_params.Jex) ) ^ ( (xdx - 1) + (ydx - 1) )
            vector *= checkerboard_stagger( model_params.Jex, xdx, ydx )
            ϕ, z = 2 * π * (-1. + 2. * rand()), -1. + 2. * rand()
            vector += model_params.initial_randomness * Spin3( cos(ϕ) * sqrt(1 - z^2), sin(ϕ) * sqrt(1 - z^2), z )
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
From an h field, calculate the ouptut spin 
"""
avg_spin( heff, β ) = unit_spin3( heff ) * LangevinFunction( abs(heff), β )

"""
Calculate the MFT spin at the site 
"""
function mft_spin_per_site( site, lattice_spins, params::ModelParameters, nearest_neighbors, one_d )
    eff_field = effective_field_per_site(site, lattice_spins, params, nearest_neighbors, one_d)
    output = avg_spin( eff_field, params.β )
    return output
end

"""
Calclate the MFT energy at a site as Eᵢ = 𝐡effⁱ ⋅ 𝐒ᵢ
"""
mft_energy_per_spin( eff_field, spin ) = -1. * eff_field ⋅ spin

"""
Sweep the lattice and compute the total energy
"""
function mft_energy_of_system( lattice_spins, params::ModelParameters, latt_params , neighbors, one_d )
    energy = 0.
    for site ∈ 1:length(lattice_spins)
        if boundary_neighbor_value != neighbors[site, 1]
            eff_field = effective_field_per_site( site, lattice_spins, params, neighbors, one_d )
            energy += mft_energy_per_spin( eff_field, lattice_spins[site] )
        end
    end
    return energy / total_bulk_sites(latt_params)
end

"""
Calculate MFT for the lattice from an initial guess
"""
function mft_lattice( lattice_spins, model_params::ModelParameters, latt_params::LatticeParameters, nearest_neighbors; iteration_scheme = nothing )
    new_spins = copy(lattice_spins)
    site_list = iteration_scheme
    if site_list === nothing
        site_list = 1:total_sites(latt_params)
    end
    for site ∈ site_list
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
        error += abs(spin)
    end
    return error / length(diff)
end