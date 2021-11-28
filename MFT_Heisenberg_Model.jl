"""
This code runs a fixed-point iteration code for 
the Heisenberg model.
"""

using Revise

include("HeisenbergModel.jl")
include("FixedPointIteration.jl")
include("PlotSpins.jl")

"""
Determine whether the model is antiferromagnetic (AFM)
or ferromagnetic (FM).
"""
function model_name( model_params::ModelParameters, latt_params::LatticeParameters ) 
    name = ""
    if model_params.Jex > 0.
        name = "AFM"
    elseif model_params.Jex < 0.
        name = "FM"
    end
    return "$(name)_Lx-$(latt_params.Lx)_Ly-$(latt_params.Ly)"
end        

function local_mft_heisenberg_main()
    figure_directory = raw"C:\Users\meese\Documents\Miscellaneous Notes\Local MFT Heisenberg Spins\Figures"

    square_L = 64
    latt_params  = LatticeParameters( square_L, square_L )
    model_params = ModelParameters( 1., 1000., 0.2 )

    nearest_neighbors = nearest_neighbor_table( latt_params )
    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    initialize_spins!(lattice_spins, latt_params, model_params)
    iteration_scheme = xy_plane_iteration_x_boundaries(latt_params)
    # iteration_scheme = nothing

    @time mft_spins, errors = FixedPointIteration( (x, y, z, w) -> mft_lattice(x, y, z, w; iteration_scheme = iteration_scheme), 
                                                  average_spin_difference, lattice_spins,
                                                  model_params, latt_params, nearest_neighbors )

    plot_spin_chain(div(latt_params.Ly, 2), latt_params, mft_spins; 
                    model_name=model_name(model_params, latt_params), save_location=figure_directory)
    plot_error_evolution( errors; 
                          model_name=model_name(model_params, latt_params), save_location=figure_directory)
    # plot_spin_arrows(latt_params, mft_spins)
    plot_spin_colormap(latt_params, mft_spins; 
                       model_name=model_name(model_params, latt_params), save_location=figure_directory)
end

@time local_mft_heisenberg_main()