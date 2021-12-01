"""
This code runs a fixed-point iteration code for 
the Heisenberg model.
"""

using Revise

include("src/HeisenbergModel.jl")
include("src/FixedPointIteration.jl")
include("src/PlotSpins.jl")
include("src/ModelNamingUtilities.jl")  

function local_mft_heisenberg_main()
    figure_directory = raw"C:\Users\meese\Documents\Miscellaneous Notes\Local MFT Heisenberg Spins\Figures"
    figure_directory = nothing

    square_L = 20
    latt_params  = LatticeParameters( square_L, square_L )
    model_params = ModelParameters( -1., 1000., 10.)

    nearest_neighbors = nearest_neighbor_table( latt_params )
    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    initialize_spins!(lattice_spins, latt_params, model_params)
    iteration_scheme = xy_plane_iteration_x_boundaries(latt_params)
    state_function = x -> mft_energy_of_system( x, model_params, latt_params, nearest_neighbors, latt_params.Ly == 1 ) 

    @time mft_spins, errors, energies = FixedPointIteration( (x, y, z, w) -> mft_lattice(x, y, z, w; iteration_scheme = iteration_scheme), 
                                                              average_spin_difference, lattice_spins,
                                                              model_params, latt_params, nearest_neighbors; 
                                                              maxiter=100000,
                                                              state_function = state_function )

    plot_spin_chain(div(latt_params.Ly, 2), latt_params, mft_spins; 
                    model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_error_evolution( errors; 
                          model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_energy_evolution( energies, errors; 
                           model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_spin_arrows(latt_params, mft_spins; chains=true)
    plot_spin_colormap(latt_params, mft_spins; 
                       model_name=model_name(model_params, latt_params), save_location=figure_directory)
end

@time local_mft_heisenberg_main()