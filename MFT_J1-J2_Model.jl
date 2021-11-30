"""
This code runs a fixed-point iteration code for 
the Heisenberg J₁-J₂ model.
"""

using Revise

include("src/Heisenberg_J1-J2_Model.jl")
include("src/FixedPointIteration.jl")
include("src/PlotSpins.jl")
include("src/ModelNamingUtilities.jl")

function local_mft_J1_J2_main()
    figure_directory = raw"C:\Users\meese\Documents\Miscellaneous Notes\Local MFT Heisenberg Spins\Figures"
    figure_directory = nothing 

    square_L = 64
    latt_params  = LatticeParameters( square_L, square_L )
    model_params = J1_J2_ModelParameters( ModelParameters(0., 1000., 0.01 ), 1. )

    nearest_neighbors  = nearest_neighbor_table( latt_params )
    Nnearest_neighbors = next_nearest_neighbor_table( latt_params )
    neigbors = ( nearest_neighbors, Nnearest_neighbors )

    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    initialize_spins!(lattice_spins, latt_params, model_params)
    iteration_scheme = xy_plane_iteration_x_boundaries(latt_params)
    # iteration_scheme = nothing

    @time mft_spins, errors = FixedPointIteration( (x, y, z, w) -> mft_lattice(x, y, z, w; iteration_scheme = iteration_scheme), 
                                                  average_spin_difference, lattice_spins,
                                                  model_params, latt_params, neigbors )

    plot_spin_chain(div(latt_params.Ly, 2), latt_params, mft_spins; 
                    model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_error_evolution( errors; 
                          model_name=model_name(model_params, latt_params), save_location=figure_directory)
    # plot_spin_arrows(latt_params, mft_spins)
    plot_spin_colormap(latt_params, mft_spins; 
                       model_name=model_name(model_params, latt_params), save_location=figure_directory)
end

@time local_mft_J1_J2_main()