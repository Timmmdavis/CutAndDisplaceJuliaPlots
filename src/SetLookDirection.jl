function SetLookDirection(scene,xmid,ymid,zmid,Prange,direction)
	#Prange 

	#Get the camera
	cam = Makie.cameracontrols(scene)
	
	if direction=="x"
		cam.eyeposition[] = (Prange.*2, ymid, 		zmid )
	elseif direction=="y"
		cam.eyeposition[] = (xmid, 		Prange.*2, 	zmid )
	elseif direction=="z"
		cam.eyeposition[] = (xmid, 		ymid, 		Prange.*2 )
	end

	cam.lookat[] = (xmid, ymid, zmid )
	update_cam!(scene, cam)
	# stop scene display from centering, which would overwrite the camera paramter we just set
	scene.center = false


	return scene

end