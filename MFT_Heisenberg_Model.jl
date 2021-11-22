"""
This code runs a fixed-point iteration code for 
the Heisenberg model.
"""

include("HeisenbergModel.jl")

latt_params  = LatticeParameters( 6, 6 )
model_params = ModelParameters( -1., 10. )

nearest_neighbors = nearest_neighbor_table( latt_params )
lattice_spins = Array{Spin3}( undef, total_sites( latt_params ) )
initialize_spins!(lattice_spins, latt_params)