"""
This file will set up the neighbor functions for the 2D 
square lattice strip. Here the system is periodic in y but
closed in x. 

Additionally, there are four boundary spins for each 
constant value of y -- two are on each side to help 
demonstrate the stripes.
"""

const num_nearest_neighbors = 4
const num_next_nearest_neighbors = 4
const boundary_neighbor_value = -2
const num_boundary_x_per_side = 2

struct LatticeParameters
    Lx::Int
    Ly::Int
end

struct Site2D
    xind::Int
    yind::Int
end

total_sites( latt_params::LatticeParameters ) = latt_params.Lx * latt_params.Ly
total_boundary_sites( latt_params::LatticeParameters ) = 2 * num_boundary_x_per_side * latt_params.Ly
total_bulk_sites( latt_params::LatticeParameters ) = total_sites(latt_params) - total_boundary_sites(latt_params)
site_index( site::Site2D, latt_params::LatticeParameters ) = (site.yind - 1) * latt_params.Lx + site.xind
function site_coords( index::Int, latt_params::LatticeParameters ) 
    xind = mod(index, latt_params.Lx)
    endtrue = xind == 0
    xind = xind * convert(Int, !endtrue) + latt_params.Lx * convert(Int, endtrue)
    return Site2D( xind, 1 + div( index - xind, latt_params.Lx ) )
end

"""
Construct the nearest neighbor table 
"""
function nearest_neighbor_table( latt_params::LatticeParameters ; num_neighbors = num_nearest_neighbors ) 
    neighbors = Array{Int}(undef, total_sites(latt_params), num_neighbors)
    for site ∈ 1:total_sites(latt_params)
        coords = site_coords( site, latt_params )
        # Populate the neighbors for the boundaries with the boundary_neighbor_value
        if coords.xind <= num_boundary_x_per_side || coords.xind > latt_params.Lx - num_boundary_x_per_side
            for nn ∈ 1:num_neighbors 
                neighbors[site, nn] = boundary_neighbor_value 
            end
        else
            # Otherwise proceed as planned. x neighbors first, then y neighbors
            neighbors[site, 1] = site_index( Site2D( coords.xind - 1, coords.yind ), latt_params )
            neighbors[site, 2] = site_index( Site2D( coords.xind + 1, coords.yind ), latt_params )
            # Cover the periodicity
            neighbor_y = (coords.yind + 1) * convert(Int, coords.yind != latt_params.Ly) + 1 * convert(Int, coords.yind == latt_params.Ly)
            neighbors[site, 3] = site_index( Site2D( coords.xind, neighbor_y ), latt_params )
            neighbor_y = (coords.yind - 1) * convert(Int, coords.yind != 1) + latt_params.Ly * convert(Int, coords.yind == 1)
            neighbors[site, 4] = site_index( Site2D( coords.xind, neighbor_y ), latt_params )
        end
    end
    return neighbors
end

"""
Construct the next-nearest neighbor table 
"""
function next_nearest_neighbor_table( latt_params::LatticeParameters ; num_neighbors = num_next_nearest_neighbors ) 
    neighbors = Array{Int}(undef, total_sites(latt_params), num_neighbors)
    for site ∈ 1:total_sites(latt_params)
        coords = site_coords( site, latt_params )
        # Populate the neighbors for the boundaries with the boundary_neighbor_value
        if coords.xind <= num_boundary_x_per_side || coords.xind > latt_params.Lx - num_boundary_x_per_side
            for nn ∈ 1:num_neighbors 
                neighbors[site, nn] = boundary_neighbor_value 
            end
        else
            # Otherwise proceed as planned. The periodicity in y must be accounted for.
            # Go backwards first. 
            neighbor_y = (coords.yind - 1) * convert(Int, coords.yind != 1) + latt_params.Ly * convert(Int, coords.yind == 1)
            neighbors[site, 1] = site_index( Site2D( coords.xind - 1, neighbor_y ), latt_params )
            neighbors[site, 2] = site_index( Site2D( coords.xind + 1, neighbor_y ), latt_params )
            # Now go forwards along y
            neighbor_y = (coords.yind + 1) * convert(Int, coords.yind != latt_params.Ly) + 1 * convert(Int, coords.yind == latt_params.Ly)
            neighbors[site, 3] = site_index( Site2D( coords.xind - 1, neighbor_y ), latt_params )
            neighbors[site, 4] = site_index( Site2D( coords.xind + 1, neighbor_y ), latt_params )
        end
    end
    return neighbors
end

neighbor = nearest_neighbor_table(LatticeParameters(5, 2))