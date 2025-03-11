function SetLookDirection(fig, xmid, ymid, zmid, Prange, direction)
    
	#= 
    # Get the Axis3 directly
    ax = fig.content[1]
    
    # Get the camera controls
    cam = ax.scene.camera
    
    # Set eye position based on direction
    if direction == "x"
        cam.eyeposition[] = [Prange * 2, ymid, zmid]
    elseif direction == "y"
        cam.eyeposition[] = [xmid, Prange * 2, zmid]
    elseif direction == "z"
        cam.eyeposition[] = [xmid, ymid, Prange * 2]
    else
        error("Invalid direction. Use 'x', 'y', or 'z'.")
    end
    
    # Set look-at point
    cam.lookat[] = [xmid, ymid, zmid]
    
    # Update camera
    update_cam!(ax.scene, cam)
    
    return fig
     =#
end