"""
This code runs a fixed-point iteration code for 
the Heisenberg J₁-J₂ model with magnetoelastic
coupling, a nematic anisotropy, and uniaxial
anisotropy.
"""

using Revise

include("src/Heisenberg_4-State_Stripes.jl")
include("src/FixedPointIteration.jl")
include("src/PlotSpins.jl")
include("src/ModelNamingUtilities.jl")

function local_mft_4_State_Stripes_main()
    figure_directory = raw"C:\Users\meese\Documents\Miscellaneous Notes\Local MFT Heisenberg Spins\Figures"
    figure_directory = nothing 

    square_L = 11
    latt_params  = LatticeParameters( square_L, square_L )
    model_params = MagElastic_Stripe_Params( J1_J2_ModelParameters( ModelParameters(0., 1000., 100.0),
                                                                    1.0 ), 0.2, 1., -0.2 )

    nearest_neighbors  = nearest_neighbor_table( latt_params )
    Nnearest_neighbors = next_nearest_neighbor_table( latt_params )
    neighbors = ( nearest_neighbors, Nnearest_neighbors )

    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    initialize_spins!(lattice_spins, latt_params, model_params)
    iteration_scheme = xy_plane_iteration_x_boundaries(latt_params)
    # iteration_scheme = nothing
    state_function = x -> mft_energy_of_system( x, model_params, latt_params, neighbors, latt_params.Ly == 1 ) 


    @time mft_spins, errors, energies = FixedPointIteration( (x, y, z, w) -> mft_lattice(x, y, z, w; iteration_scheme = iteration_scheme), 
                                                  average_spin_difference, lattice_spins,
                                                  model_params, latt_params, neighbors; 
                                                  maxiter = 10000,
                                                  state_function=state_function )

    plot_spin_chain(div(latt_params.Ly, 2), latt_params, mft_spins; 
                    model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_error_evolution( errors; 
                          model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_energy_evolution( energies, errors; 
                           model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_spin_arrows(latt_params, mft_spins; chains=true)
    plot_spin_colormap(latt_params, mft_spins; 
                       model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_function_of_x( latt_params.Lx, x -> nematicity.(x, model_params.ε, latt_params.Lx ), "\$\\mathrm{Strain}\$ \$\\varepsilon(x)\$" )
end

@time local_mft_4_State_Stripes_main()