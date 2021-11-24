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
    
    fig, ax = PyPlot.subplots(3,1, figsize=(6,6), sharex=true, sharey=true)
    ax[1].plot( mft_S[:,1], marker="o", mec="k" )
    ax[1].set_ylabel(L"$\left\langle S^x(x) \right\rangle$")
    
    ax[2].plot( mft_S[:,2], marker="o", mec="k" )
    ax[2].set_ylabel(L"$\left\langle S^y(x) \right\rangle$")
    
    ax[3].plot( mft_S[:,3], marker="o", mec="k" )
    ax[3].set_ylabel(L"$\left\langle S^z(x) \right\rangle$")

    ax[3].set_xlabel(L"Site along chain $x$")
    fig.tight_layout()
    PyPlot.show()
end

function local_mft_heisenberg_main()
    latt_params  = LatticeParameters( 7, 1 )
    model_params = ModelParameters( -1., 1000. )

    nearest_neighbors = nearest_neighbor_table( latt_params )
    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    @time initialize_spins!(lattice_spins, latt_params)

    lattice_spins

    @time mft_spins = FixedPointIteration(mft_lattice, average_spin_difference, lattice_spins,
                                        model_params, latt_params, nearest_neighbors )

    plot_spin_chain(1, latt_params, mft_spins)
end

@time local_mft_heisenberg_main()