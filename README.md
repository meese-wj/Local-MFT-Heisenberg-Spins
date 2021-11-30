# Local Mean Field Theory for a Heisenberg Chain

This repo contains a toy code that uses fixed-point iteration to obtain a mean-field solution to a 2d Heisenberg model.
The mean-field value of each spin is obtained per site using this approach on a 2d square lattice which is periodic in
the y-direction but has fixed boundaries in the x-direction.

## Important `Julia` scripts
1. The `MFT_Heisenberg_Model.jl` file runs a nearest-neighbor system of size Lx by Ly. The x-direction has fixed boundary 
conditions which need to be specified, and the y-direction has periodic boundary conditions.

1. The `MFT_J1-J2_Model.jl` file runs a system with nearest and next-nearest neighbor bonds, on a lattice that looks like 
that in 1.