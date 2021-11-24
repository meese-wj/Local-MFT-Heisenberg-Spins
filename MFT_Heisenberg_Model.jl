"""
This code runs a fixed-point iteration code for 
the Heisenberg model.
"""

using Revise
using PyPlot

include("HeisenbergModel.jl")
include("FixedPointIteration.jl")

function local_mft_heisenberg_main()
    latt_params  = LatticeParameters( 200, 1 )
    model_params = ModelParameters( -1., 10. )

    nearest_neighbors = nearest_neighbor_table( latt_params )
    lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
    @time initialize_spins!(lattice_spins, latt_params)

    lattice_spins

    @time mft_spins = FixedPointIteration(mft_lattice, average_spin_difference, lattice_spins,
                                        model_params, latt_params, nearest_neighbors )

    mft_S₃ = zeros( latt_params.Lx )
    for xdx ∈ 1:latt_params.Lx
        mft_S₃[xdx] = mft_spins[ site_index( Site2D(xdx, 1), latt_params ) ].S₃
    end
    PyPlot.plot( mft_S₃ )
    PyPlot.show()
end

@time local_mft_heisenberg_main()