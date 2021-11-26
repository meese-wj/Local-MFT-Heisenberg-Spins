"""
This file defines a bunch of different 
plotting routines for Heisenberg spins.
"""
using Revise
using PyPlot

PyPlot.rc("xtick", direction="in")
PyPlot.rc("ytick", direction="in")

include("LatticeSetup.jl")
include("HeisenbergSpins.jl")

"""
Extract the spin components and put them in
an (N x 3) array.
"""
function spins_to_array( latt_params, mft_spins )
    mft_S = zeros( latt_params.Ly, latt_params.Lx, 3 )
    for ydx ∈ 1:latt_params.Ly, xdx ∈ 1:latt_params.Lx
        mft_S[ydx, xdx, 1] = mft_spins[ site_index( Site2D(xdx, ydx), latt_params ) ].S₁
        mft_S[ydx, xdx, 2] = mft_spins[ site_index( Site2D(xdx, ydx), latt_params ) ].S₂
        mft_S[ydx, xdx, 3] = mft_spins[ site_index( Site2D(xdx, ydx), latt_params ) ].S₃
    end
    return mft_S
end

"""
Plot a shaded region for the bounday spins. This is 
essentially a subroutine to clean up code.
"""
function plot_boundary_spins( ax, num_boundary, latt_params )
    ymin, ymax = ax[3].get_ylim()
    for component ∈ 1:3
        left_boundary = LinRange(1, num_boundary, 10)
        right_boundary = LinRange(latt_params.Lx + 1 - num_boundary, latt_params.Lx, 10)
        ax[component].fill_between( left_boundary, ymin .+  0 .* left_boundary , ymax .+ 0 .* left_boundary, color = "orange", alpha=0.3 )
        ax[component].fill_between( right_boundary, ymin .+  0 .* right_boundary , ymax .+ 0 .* right_boundary, color = "orange", alpha=0.3 )
    end
    return ax
end

"""
Plot a single spin chain with fixed yindex.
"""
function plot_spin_chain( yindex, latt_params, mft_spins )
    mft_S = spins_to_array(latt_params, mft_spins)
    
    xvalues = LinRange(1, latt_params.Lx, latt_params.Lx)
    marker = "o"
    if latt_params.Lx > 50
        marker = "None"
    end

    fig, ax = PyPlot.subplots(3,1, figsize=(6,6), sharex=true, sharey=true)
    ax[1].plot( xvalues, mft_S[yindex, :, 1], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[1].set_ylabel(L"$\left\langle S^x(x) \right\rangle$")
    
    ax[2].plot( xvalues, mft_S[yindex, :, 2], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[2].set_ylabel(L"$\left\langle S^y(x) \right\rangle$")
    
    ax[3].plot( xvalues, mft_S[yindex, :, 3], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[3].set_ylabel(L"$\left\langle S^z(x) \right\rangle$")

    ax = plot_boundary_spins(ax, 2, latt_params)

    ax[3].set_xlabel(L"Site along chain $x$")
    ax[3].set_xlim(1, latt_params.Lx)
    ax[3].set_ylim(ymin, ymax)
    fig.tight_layout()
    PyPlot.show()
end

"""
Plot the fixed-point iteration algorithm error
as a function of the iteration.
"""
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