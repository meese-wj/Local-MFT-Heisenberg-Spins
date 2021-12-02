"""
This file defines a bunch of different 
plotting routines for Heisenberg spins.
"""

using Revise
using PyPlot
using PyCall
@pyimport matplotlib.colors as mpl_colors

PyPlot.rc("xtick", direction="in")
PyPlot.rc("ytick", direction="in")
PyPlot.rc("font", size=12)        # Set this as the default size
# PyPlot.rc("text",  usetex=true)   # This is slow upon startup due to TeX

include("LatticeSetup.jl")
include("HeisenbergSpins.jl")

""" 
Write a wrapper around the PyPlot 
savefig function that also creates the
file directory housing the figure output.

If save_location is outside of the runtime
directory, then a full path should be used.
"""
function figure_save_wrapper(figure, figure_name, save_location, extension; savefig_kwargs...)
    if save_location === nothing
        return 
    end
    save_location = mkpath(save_location)
    ext = extension
    if ext[1] == '.'
        ext = extension[2:end]
    end
    fig_file = joinpath(save_location, "$figure_name.$ext")
    figure.savefig(fig_file, bbox_inches="tight", savefig_kwargs...)
    return
end

"""
Use this to add model names to the figures.
"""
function prepend_model_name( model_name::String, file_string::String )
    output = file_string
    if length(model_name) > 0
        output = "$(model_name)_$file_string"
    end
    return output
end

"""
Plot the fixed-point iteration algorithm error
as a function of the iteration.
"""
function plot_error_evolution( all_errors; 
                               model_name="", save_location=nothing, extension=".pdf" )
    first_flag = findfirst(x -> x == error_notice_flag, all_errors )
    if first_flag === nothing
        first_flag = length(all_errors)
    else
        first_flag -= 1
    end
    fig = PyPlot.figure()
    PyPlot.loglog( LinRange( 1, first_flag-1, first_flag ), all_errors[begin : first_flag], lw=3 )
    PyPlot.xlabel(L"$\mathrm{Fixed-point\, iteration}$")
    PyPlot.ylabel(L"$\mathrm{Error\, per\, site}$")
    PyPlot.grid(which="major")
    PyPlot.tight_layout()

    figure_save_wrapper( fig, prepend_model_name(model_name, "error_per_spin"), 
                         save_location, extension )
    PyPlot.show()
end

"""
Plot the fixed-point iteration algorithm error
as a function of the iteration.
"""
function plot_energy_evolution( all_energies, all_errors; 
                                model_name="", save_location=nothing, extension=".pdf" )
    first_flag = findfirst(x -> x == error_notice_flag, all_errors )
    if first_flag === nothing
        first_flag = length(all_errors)
    else
        first_flag -= 1
    end
    fig = PyPlot.figure()
    PyPlot.semilogx( LinRange( 1, first_flag-1, first_flag ), all_energies[begin : first_flag], lw=3 )
    PyPlot.xlabel(L"$\mathrm{Fixed-point\, iteration}$")
    PyPlot.ylabel(L"$\mathrm{Mean\, field\, energy}$")
    PyPlot.grid(which="major")
    PyPlot.tight_layout()

    figure_save_wrapper( fig, prepend_model_name(model_name, "error_per_spin"), 
                         save_location, extension )
    PyPlot.show()
end

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
    ymin, ymax = 0., 0.
    for axis in ax
        ymin, ymax = axis.get_ylim()
    end
    for component ∈ 1:3
        left_boundary = LinRange(1, num_boundary, 10)
        right_boundary = LinRange(latt_params.Lx + 1 - num_boundary, latt_params.Lx, 10)
        ax[component].fill_between( left_boundary, ymin .+  0 .* left_boundary , ymax .+ 0 .* left_boundary, color = "orange", alpha=0.3 )
        ax[component].fill_between( right_boundary, ymin .+  0 .* right_boundary , ymax .+ 0 .* right_boundary, color = "orange", alpha=0.3 )
    end
    for axis in ax
        axis.set_ylim(ymin, ymax)
    end
    return ax
end

"""
Plot a single spin chain with fixed yindex.
"""
function plot_spin_chain( yindex, latt_params, mft_spins; 
                          model_name="", save_location=nothing, extension=".pdf" )
    mft_S = spins_to_array(latt_params, mft_spins)
    
    xvalues = LinRange(1, latt_params.Lx, latt_params.Lx)
    marker = "o"
    if latt_params.Lx > 50
        marker = "None"
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    default_fontsize = rcParams["font.size"]
    rcParams["font.size"] = 16
    
    fig, ax = PyPlot.subplots(3,1, figsize=(6,6), sharex=true, sharey=true)
    ax[1].set_ylim(-1., 1.)
    ax[1].plot( xvalues, mft_S[yindex, :, 1], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[1].set_ylabel("\$\\left\\langle S^x(x, $yindex) \\right\\rangle\$")
    
    ax[2].plot( xvalues, mft_S[yindex, :, 2], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[2].set_ylabel("\$\\left\\langle S^y(x, $yindex) \\right\\rangle\$")
    
    ax[3].plot( xvalues, mft_S[yindex, :, 3], marker=marker, mec="k", clip_on=false, zorder=20 )
    ax[3].set_ylabel("\$\\left\\langle S^z(x, $yindex) \\right\\rangle\$")

    ax = plot_boundary_spins(ax, num_boundary_x_per_side, latt_params)

    ax[3].set_xlabel(L"$x$ $\mathrm{site\, along\, chain}$")
    ax[3].set_xlim(1, latt_params.Lx)
    fig.tight_layout()

    figure_save_wrapper( fig, prepend_model_name(model_name, "LMFT_chain_projections"), 
                         save_location, extension )
    PyPlot.show()
    rcParams["font.size"] = default_fontsize
end

"""
Plot the spins as arrows in 3D
"""
function plot_spin_arrows(latt_params, mft_spins; chains=false,
                          model_name="", save_location=nothing, extension=".pdf")
    z_plane_coord = 0.
    xyz_coords = zeros( total_sites(latt_params), 3 )
    for site ∈ 1:total_sites(latt_params)
        coord = site_coords( site, latt_params )
        xyz_coords[site, 1] = coord.xind
        xyz_coords[site, 2] = coord.yind
        xyz_coords[site, 3] = z_plane_coord
    end

    fig = PyPlot.figure()
    ax = fig.add_subplot(projection="3d")
    arrow_length = latt_params.Lx / 10.
    ax.set_box_aspect((1,1,1))
    if !chains
        for site ∈ 1:total_sites(latt_params)
            ax.quiver( xyz_coords[site, 1], xyz_coords[site, 2], xyz_coords[site, 3],
                    mft_spins[site].S₁, mft_spins[site].S₂, mft_spins[site].S₃,
                    arrow_length_ratio=0.15, length=arrow_length )
        end
    else
        yindices = [ div(latt_params.Lx, 4) div(latt_params.Lx, 2) div(3 * latt_params.Lx, 4) ]
        colors = ["red" "green" "blue"]
        for (cdx, ydx) ∈ enumerate(yindices), xdx ∈ 1:latt_params.Lx
            site = site_index( Site2D(xdx, ydx), latt_params )
            ax.quiver3D( xyz_coords[site, 1], xyz_coords[site, 2], xyz_coords[site, 3],
                    mft_spins[site].S₁, mft_spins[site].S₂, mft_spins[site].S₃,
                    arrow_length_ratio=0.15, length=arrow_length, color=colors[cdx] )
            
        end
    end
    ax.set_xlabel(L"$x$")
    ax.set_ylabel(L"$y$")
    ax.set_zlabel(L"$z$")
    ax.set_xlim(0, latt_params.Lx)
    ax.set_ylim(0, latt_params.Lx)
    ax.set_zlim(-div(latt_params.Lx, 2), div(latt_params.Lx, 2))
    ax.grid(false)

    figure_save_wrapper( fig, prepend_model_name(model_name, "LMFT_xy_plane_arrows"), 
                         save_location, extension )
    PyPlot.show()    
end

# Define a colormap for use in plotting the spins
my_gradient = mpl_colors.LinearSegmentedColormap.from_list("my_gradient", (
                                                            # Edit this gradient at https://eltos.github.io/gradient/#4C71FF-0025B3-000000-C7030D-FC4A53
                                                            (0.000, (0.298, 0.443, 1.000)),
                                                            (0.250, (0.000, 0.145, 0.702)),
                                                            (0.500, (0.000, 0.000, 0.000)),
                                                            (0.750, (0.780, 0.012, 0.051)),
                                                            (1.000, (0.988, 0.290, 0.325))))


blue_orange_gradient = mpl_colors.LinearSegmentedColormap.from_list("blue_orange_gradient", (
                                                                    # Edit this gradient at https://eltos.github.io/gradient/#4C71FF-0025B3-000000-B74A01-FF902C
                                                                    (0.000, (0.298, 0.443, 1.000)),
                                                                    (0.250, (0.000, 0.145, 0.702)),
                                                                    (0.500, (0.000, 0.000, 0.000)),
                                                                    (0.750, (0.780, 0.431, 0.012)),
                                                                    (1.000, (0.988, 0.659, 0.290))))

"""
Plot spin colormap
"""
function plot_spin_colormap(latt_params, mft_spins; model_name="", save_location=nothing, extension=".pdf")
    Sx_values = zeros(latt_params.Lx, latt_params.Ly)
    Sy_values = zeros(latt_params.Lx, latt_params.Ly)
    Sz_values = zeros(latt_params.Lx, latt_params.Ly)
    for site ∈ 1:total_sites(latt_params)
        coords = site_coords(site, latt_params)
        Sx_values[coords.xind, coords.yind] = mft_spins[site].S₁
        Sy_values[coords.xind, coords.yind] = mft_spins[site].S₂
        Sz_values[coords.xind, coords.yind] = mft_spins[site].S₃
    end

    cbar_min = min( min(minimum(Sx_values), minimum(Sy_values)), minimum(Sz_values) )
    cbar_max = max( max(maximum(Sx_values), maximum(Sy_values)), maximum(Sz_values) )
    the_cmap = blue_orange_gradient

    fig, axs = PyPlot.subplots(1, 3, sharex=true, sharey=true)
    axs[1].imshow( Sx_values', origin="lower", vmin = cbar_min, vmax = cbar_max, cmap=the_cmap, extent=[1, latt_params.Lx, 1, latt_params.Ly] )
    axs[1].set_title(L"$\left\langle S^x(x,y)\right\rangle$")
    axs[1].set_ylabel(L"$y$")
    axs[1].set_xlabel(L"$x$")
    axs[2].imshow( Sy_values', origin="lower", vmin = cbar_min, vmax = cbar_max, cmap=the_cmap, extent=[1, latt_params.Lx, 1, latt_params.Ly] )
    axs[2].set_xlabel(L"$x$")
    axs[2].set_title(L"$\left\langle S^y(x,y)\right\rangle$")
    sz_im = axs[3].imshow( Sz_values', origin="lower", vmin = cbar_min, vmax = cbar_max, cmap=the_cmap, extent=[1, latt_params.Lx, 1, latt_params.Ly] )
    axs[3].set_xlabel(L"$x$")
    axs[3].set_title(L"$\left\langle S^z(x,y)\right\rangle$")

    for comp ∈ 1:3
        axs[comp].tick_params(which="both", bottom=false, left=false)
    end

    fig.tight_layout()
    sx_bbox = axs[1].get_position().bounds  # (x, y, width, height)
    sz_bbox = axs[3].get_position().bounds 
    
    # Now add the colorbar
    fig.subplots_adjust(bottom=0.25)  # Move the bottom of the subplots up by 20%
    cbar_width = sz_bbox[1] + sz_bbox[3] - sx_bbox[1]
    cbar_axis = fig.add_axes([sx_bbox[1], 0.15, cbar_width, 0.05])
    cbar = fig.colorbar(sz_im, cax=cbar_axis, orientation="horizontal")
    cbar.set_label(L"\rm Spin\, Projection", loc="center")

    # Crop the figure
    fig_height = fig.get_figheight()
    new_figheight = (sz_bbox[2] + sz_bbox[4]) * fig_height
    fig.set_figheight(new_figheight)

    figure_save_wrapper( fig, prepend_model_name(model_name, "LMFT_xy_plane_projections"), 
                         save_location, extension )

    PyPlot.show()
end

function plot_function_of_x( Lx, func::Function, ylabel::String )
    xvalues = LinRange(1, Lx, Lx)
    yvalues = func( xvalues )

    fig, ax = PyPlot.subplots(1,1)
    ax.plot(xvalues, yvalues)
    ax.set_xlim(1, Lx)
    ax.set_xlabel(L"$x$ $\mathrm{site\, along\, chain}$")
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    return
end

