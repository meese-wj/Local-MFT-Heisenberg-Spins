"""
This code runs a fixed-point iteration code for 
the Heisenberg model.
"""

using Revise
using PyPlot

include("HeisenbergModel.jl")
include("FixedPointIteration.jl")

function plot_spin_chain( yindex, latt_params, mft_spins )
    mft_S = zeros( latt_params.Lx, 3 )
    for xdx ∈ 1:latt_params.Lx
        mft_S[xdx, 1] = mft_spins[ site_index( Site2D(xdx, yindex), latt_params ) ].S₁
        mft_S[xdx, 2] = mft_spins[ site_index( Site2D(xdx, yindex), latt_params ) ].S₂
        mft_S[xdx, 3] = mft_spins[ site_index( Site2D(xdx, yindex), latt_params ) ].S₃
    end
    
    xvalues = LinRange(1, latt_params.Lx, latt_params.Lx)

    fig, ax = PyPlot.subplots(3,1, figsize=(6,6), sharex=true, sharey=true)
    ax[1].plot( xvalues, mft_S[:,1], marker="o", mec="k", clip_on=false, zorder=20 )
    ax[1].set_ylabel(L"$\left\langle S^x(x) \right\rangle$")
    
    ax[2].plot( xvalues, mft_S[:,2], marker="o", mec="k", clip_on=false, zorder=20 )
    ax[2].set_ylabel(L"$\left\langle S^y(x) \right\rangle$")
    
    ax[3].plot( xvalues, mft_S[:,3], marker="o", mec="k", clip_on=false, zorder=20 )
    ax[3].set_ylabel(L"$\left\langle S^z(x) \right\rangle$")

    ax[3].set_xlabel(L"Site along chain $x$")
    ax[3].set_xlim(1, latt_params.Lx)
    fig.tight_layout()
    PyPlot.show()
end

function plot_error_evolution( all_errors )
    first_flag = findfirst(x -> x == error_notice_flag, all_errors )
    if first_flag === nothing
        first_flag = length(all_errors)
    end
    PyPlot.figure()
    PyPlot.loglog( LinRange( 1, first_flag-1, first_flag ), all_errors[begin : first_flag], lw=3 )
    PyPlot.xlabel("Fixed-point iteration Step")
    PyPlot.ylabel("Step-wise error per site")
    PyPlot.grid(which="both")
    PyPlot.tight_layout()
    PyPlot.show()
end

function local_mft_heisenberg_main()
    latt_params  = LatticeParameters( 10, 3 )
    model_params = ModelParameters( -1., 1000., 0.01 )

    nearest_neighbors = nearest_neighbor_table( latt_params )
    display(nearest_neighbors)
    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    @time initialize_spins!(lattice_spins, latt_params, model_params)

    lattice_spins

    @time mft_spins, errors = FixedPointIteration(mft_lattice, average_spin_difference, lattice_spins,
                                                  model_params, latt_params, nearest_neighbors )

    plot_spin_chain(1, latt_params, mft_spins)
    plot_error_evolution( errors )
end

@time local_mft_heisenberg_main()