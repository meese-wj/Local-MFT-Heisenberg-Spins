"""
This file defines a bunch of different 
plotting routines for Heisenberg spins.
"""
using Revise
using PyPlot

PyPlot.rc("xtick", direction="in")
PyPlot.rc("ytick", direction="in")

function plot_spin_chain( yindex, latt_params, mft_spins )
    mft_S = zeros( latt_params.Lx, 3 )
    for xdx ∈ 1:latt_params.Lx
        mft_S[xdx, 1] = mft_spins[ site_index( Site2D(xdx, yindex), latt_params ) ].S₁
        mft_S[xdx, 2] = mft_spins[ site_index( Site2D(xdx, yindex), latt_params ) ].S₂
        mft_S[xdx, 3] = mft_spins[ site_index( Site2D(xdx, yindex), latt_params ) ].S₃
    end
    
    xvalues = LinRange(1, latt_params.Lx, latt_params.Lx)
    marker = "o"
    if latt_params.Lx > 50
        marker = "None"
    end

    fig, ax = PyPlot.subplots(3,1, figsize=(6,6), sharex=true, sharey=true)
    ax[1].plot( xvalues, mft_S[:,1], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[1].set_ylabel(L"$\left\langle S^x(x) \right\rangle$")
    
    ax[2].plot( xvalues, mft_S[:,2], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[2].set_ylabel(L"$\left\langle S^y(x) \right\rangle$")
    
    ax[3].plot( xvalues, mft_S[:,3], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[3].set_ylabel(L"$\left\langle S^z(x) \right\rangle$")

    ymin, ymax = ax[3].get_ylim()
    for component ∈ 1:3
        left_boundary = LinRange(1, 2, 10)
        right_boundary = LinRange(latt_params.Lx - 1, latt_params.Lx, 10)
        ax[component].fill_between( left_boundary, ymin .+  0 .* left_boundary , ymax .+ 0 .* left_boundary, color = "orange", alpha=0.3 )
        ax[component].fill_between( right_boundary, ymin .+  0 .* right_boundary , ymax .+ 0 .* right_boundary, color = "orange", alpha=0.3 )
    end

    ax[3].set_xlabel(L"Site along chain $x$")
    ax[3].set_xlim(1, latt_params.Lx)
    ax[3].set_ylim(ymin, ymax)
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