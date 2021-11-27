"""
This code runs a fixed-point iteration code for 
the Heisenberg model.
"""

using Revise

include("HeisenbergModel.jl")
include("FixedPointIteration.jl")
include("PlotSpins.jl")

function local_mft_heisenberg_main()
    square_L = 6
    latt_params  = LatticeParameters( square_L, 3 )
    model_params = ModelParameters( 1., 1000., 0.15 )

    nearest_neighbors = nearest_neighbor_table( latt_params )
    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    @time initialize_spins!(lattice_spins, latt_params, model_params)
    iteration_scheme = xy_plane_iteration_x_boundaries(latt_params)
    display(iteration_scheme)

    @time mft_spins, errors = FixedPointIteration( x -> mft_lattice(x, model_params, latt_params, nearest_neighbors; iteration_scheme = iteration_scheme), 
                                                  average_spin_difference, lattice_spins,
                                                  model_params, latt_params, nearest_neighbors )

    plot_spin_chain(1, latt_params, mft_spins)
    plot_error_evolution( errors )
    # plot_spin_arrows(latt_params, mft_spins)
    plot_spin_colormap(latt_params, mft_spins)
end

@time local_mft_heisenberg_main()