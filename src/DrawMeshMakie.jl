function DrawMeshMakie(P1, P2, P3)
    # Draws a triangulated mesh with Makie (08/4/19) showing only the boundaries
    # P1 P2 P3 - Tri points in a list. Each row (of all three) representing a tri

    # Preparing mesh data
    x = vec([P1[:,1] P2[:,1] P3[:,1]])
    y = vec([P1[:,2] P2[:,2] P3[:,2]])
    z = vec([P1[:,3] P2[:,3] P3[:,3]])

    # List in sequential order (each 3 cols in list (row vec) is a tri)
    n = size(P1, 1)

    # Create figure
    fig = Makie.Figure()

    # Create 3D axis
    ax = Makie.Axis3(fig[1, 1], 
                aspect = :data, 
                xlabel = "X", 
                ylabel = "Y", 
                zlabel = "Z")

    # Add a wireframe for all triangles - use custom wireframe implementation to avoid dependency issues
    # Extract all edges
    edges = Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()

    for i in 1:n
        # Create three edges for each triangle
        edge1 = ([P1[i,1], P2[i,1]], [P1[i,2], P2[i,2]], [P1[i,3], P2[i,3]])
        edge2 = ([P2[i,1], P3[i,1]], [P2[i,2], P3[i,2]], [P2[i,3], P3[i,3]])
        edge3 = ([P3[i,1], P1[i,1]], [P3[i,2], P1[i,2]], [P3[i,3], P1[i,3]])

        push!(edges, edge1, edge2, edge3)
    end

    # Draw all edges manually
    for (ex, ey, ez) in edges
        Makie.lines!(ax, ex, ey, ez, color = :black, linewidth = 1.5)
    end

    # Optional: Create a semi-transparent mesh for the surface
    # For each triangle, create a separate mesh with transparency
    for i in 1:n
        # Get the three vertices of this triangle
        triangle_x = [P1[i,1], P2[i,1], P3[i,1]]
        triangle_y = [P1[i,2], P2[i,2], P3[i,2]]
        triangle_z = [P1[i,3], P2[i,3], P3[i,3]]

        # Draw triangle as a mesh with transparency
        Makie.mesh!(ax, triangle_x, triangle_y, triangle_z, 
                  color = Makie.RGBAf(0, 0, 1, 0.2),  # Blue with high transparency
                  shading = Makie.NoShading)
    end

    # Set axis limits to ensure all data is visible
    x_min, x_max = minimum(x), maximum(x)
    y_min, y_max = minimum(y), maximum(y)
    z_min, z_max = minimum(z), maximum(z)

    # Add a small buffer around the data
    buffer = 0.05
    x_range = max(x_max - x_min, 0.001)  # Ensure non-zero range
    y_range = max(y_max - y_min, 0.001)  # Ensure non-zero range
    z_range = max(z_max - z_min, 0.001)  # Ensure non-zero range

    # Set axis limits, ensuring they're not identical
    Makie.xlims!(ax, x_min - buffer * x_range, x_max + buffer * x_range)
    Makie.ylims!(ax, y_min - buffer * y_range, y_max + buffer * y_range)

    # For z-axis, make sure we have a minimum range even if the data is flat
    if z_min ≈ z_max
        # If flat in z-direction, create an artificial range
        z_center = z_min
        artificial_range = max(0.001, 0.1 * (x_range + y_range) / 2)
        Makie.zlims!(ax, z_center - artificial_range / 2, z_center + artificial_range / 2)
    else
        Makie.zlims!(ax, z_min - buffer * z_range, z_max + buffer * z_range)
    end

    return fig
end