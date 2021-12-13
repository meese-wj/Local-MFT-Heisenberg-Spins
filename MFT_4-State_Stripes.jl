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

function local_mft_4_State_Stripes_main(; λ=0.25, γ2λ=-0.5, square_Lx=8, spin_plots=true)
    figure_directory = raw"C:\Users\meese\Documents\Miscellaneous Notes\Local MFT Heisenberg Spins\Figures"
    figure_directory = nothing 

    square_L = square_Lx
    latt_params  = LatticeParameters( square_L, square_L )
    model_params = MagElastic_Stripe_Params( J1_J2_ModelParameters( ModelParameters(0.3, 1000., 100.0),
                                                                    1.0 ), λ, 0.15, γ2λ * λ )

    nearest_neighbors  = nearest_neighbor_table( latt_params )
    Nnearest_neighbors = next_nearest_neighbor_table( latt_params )
    neighbors = ( nearest_neighbors, Nnearest_neighbors )

    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    initialize_spins!(lattice_spins, latt_params, model_params)
    iteration_scheme = xy_plane_iteration_x_boundaries(latt_params)
    # iteration_scheme = nothing
    state_function = x -> mft_energy_of_system( x, model_params, latt_params, neighbors, latt_params.Ly == 1 ) 


    mft_spins, errors, energies = FixedPointIteration( (x, y, z, w) -> mft_lattice(x, y, z, w; iteration_scheme = iteration_scheme), 
                                                  average_spin_difference, lattice_spins,
                                                  model_params, latt_params, neighbors; 
                                                  maxiter = 10000,
                                                  state_function=state_function )

    if spin_plots
        plot_spin_chain(div(latt_params.Ly, 2), latt_params, mft_spins; 
                        model_name=model_name(model_params, latt_params), save_location=figure_directory)
        plot_error_evolution( errors; 
                            model_name=model_name(model_params, latt_params), save_location=figure_directory)
        plot_energy_evolution( energies, errors; 
                            model_name=model_name(model_params, latt_params), save_location=figure_directory)
        plot_spin_arrows(latt_params, mft_spins; chains=true)
        plot_spin_colormap(latt_params, mft_spins; 
                        model_name=model_name(model_params, latt_params), save_location=figure_directory)
        plot_function_of_site( latt_params.Lx, x -> nematicity.(x, model_params.ε, latt_params.Lx ), "\$\\mathrm{Strain}\$ \$\\varepsilon(x)\$" )
        # return
    end
    
    first_flag = findfirst(x -> x == error_notice_flag, errors )
    if first_flag === nothing
        first_flag = length(errors)
    else
        first_flag -= 1
    end
    return energies[first_flag]
end

function final_energies_with_system_size( γ2λ_values, λ, final_Lx )
    fig, ax = PyPlot.subplots(1,1)
    
    Lx_values = (7:final_Lx)
    final_energies = zeros( length(Lx_values), length(γ2λ_values) )

    for (γdx, γ) ∈ enumerate(γ2λ_values)
        for (idx, Ldx) ∈ enumerate(Lx_values)
            final_energies[idx, γdx] = local_mft_4_State_Stripes_main(; γ2λ=γ, square_Lx=Ldx, spin_plots=false)
        end
        
        line, = ax.plot( Lx_values, final_energies[:, γdx], label="\$\\gamma = $γ \\,\\lambda\$" )
        # ax.plot( Lx_values[ isodd.(Lx_values) ], final_energies[ isodd.(Lx_values), γdx ], color=line.get_color(), 
        #          label="\$\\gamma = $γ \\,\\lambda\$ Odd \$L_x\$", marker="o", mec="k", ls="None" )
        # ax.plot( Lx_values[ iseven.(Lx_values) ], final_energies[ iseven.(Lx_values), γdx ], 
        #          label="\$\\gamma = $γ \\,\\lambda\$ Even \$L_x\$", marker="o", mec="k", ls="None" )
    end

    for odd_idx ∈ Lx_values[isodd.(Lx_values)]
        print(odd_idx, " ")
        ax.axvline(odd_idx, lw=1, ls="dashed", color="gray", zorder=0)
    end

    ax.legend(framealpha=1.0)
    ax.set_xlabel(L"$L_x$")
    ax.set_ylabel(L"$E_{\mathrm{MF}}$ $\mathrm{per\, bulk\, spin}$ $(J_2)$")
    ax.set_title("\$ \\lambda = $λ \\, J_2 \$")
    PyPlot.show()
end

@time local_mft_4_State_Stripes_main(; square_Lx=32, spin_plots=true)
# @time final_energies_with_system_size( [-2.5 -2. -1.5 -1.0 -0.5 -0.25 -0.125], 0.2, 20)