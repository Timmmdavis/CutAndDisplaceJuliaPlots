function DrawMeshMakieFilledFaces(Val, ValName, P1, P2, P3, ColorMap)
    # Draws a triangulated mesh with Makie with a set value over the meshes faces
    # If you dont have a value set Val to: Val=ones(size(P1[:,1]))
    # And ValName to []; 
    # Val - value you want to Draw
    # ValName - Name of the value
    # P1 P2 P3 - Tri points in a list. Each row (of all three) representing a tri
    # ColorMap - The colormap you want to use
    
    # Finding respective places in the to map the colours too
    MaxVal = maximum(Val)
    MinVal = minimum(Val)
    cmapsz = size(ColorMap, 1)
    
    # Place values in 0-1 range and map to colormap indices
    Val = Val .- MinVal
    if maximum(Val) != 0  # Avoid division by zero
        Val = Val ./ maximum(Val)
    end
    Val = (Val .* (cmapsz-1)) .+ 1
    Val = round.(Val) 
    Val = convert(Array{Int64,1}, Val)
    
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
    
    # For each triangle, create a separate mesh
    for i in 1:n
        # Get the three vertices of this triangle
        triangle_x = [P1[i,1], P2[i,1], P3[i,1]]
        triangle_y = [P1[i,2], P2[i,2], P3[i,2]]
        triangle_z = [P1[i,3], P2[i,3], P3[i,3]]
        
        # Ensure index is within bounds of the colormap
        color_idx = min(max(Int(Val[i]), 1), cmapsz)
        
        # Get the color for this triangle from the colormap
        tri_color = Makie.RGB(
            ColorMap[color_idx, 1],
            ColorMap[color_idx, 2],
            ColorMap[color_idx, 3]
        )
        
        # Draw triangle as a mesh
        Makie.mesh!(ax, triangle_x, triangle_y, triangle_z, color = tri_color, shading = Makie.NoShading)
    end
    
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
        Makie.lines!(ax, ex, ey, ez, color = :black, linewidth = 1, transparency = true)
    end
    
    # Create a simple colorbar using a more direct approach
    # Create values along x axis for colorbar
    x_range = LinRange(MinVal, MaxVal, cmapsz)
    
    # Create a new axis for the colorbar
    colorbar_ax = Makie.Axis(fig[2, 1], 
                           ylabel = ValName, 
                           xlabel = "",
                           aspect = Makie.AxisAspect(10))  # Changed to AxisAspect which takes a ratio
    
    # Create a series of small colored rectangles
    for i in 1:(cmapsz-1)
        color = Makie.RGB(ColorMap[i, 1], ColorMap[i, 2], ColorMap[i, 3])
        x1, x2 = x_range[i], x_range[i+1]
        Makie.poly!(colorbar_ax, 
                  [Makie.Point2f(x1, 0), Makie.Point2f(x2, 0), 
                   Makie.Point2f(x2, 1), Makie.Point2f(x1, 1)], 
                  color = color)
    end
    
    # Set the limits of the colorbar axis
    Makie.xlims!(colorbar_ax, MinVal, MaxVal)
    Makie.ylims!(colorbar_ax, 0, 1)
    
    # Hide y-axis ticks since they're not meaningful
    Makie.hideydecorations!(colorbar_ax, grid = false)
    
    # Adjust layout
    Makie.rowsize!(fig.layout, 1, Makie.Relative(0.8))
    Makie.rowsize!(fig.layout, 2, Makie.Relative(0.2))
    
    return fig
end