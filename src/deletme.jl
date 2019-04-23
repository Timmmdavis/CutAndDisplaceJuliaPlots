#Some example plotting (Can do during debug)
scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(Dn,"Dn",P1,P2,P3,cmap2)
using Makie
#Plots to current scene if a scene is not set by default
scatter!(Px, Py, Pz, markersize = 50)

#Show normals - scale length scale by max sep of tris
arrows!(MidPoint[:,1],MidPoint[:,2],MidPoint[:,3],
	FaceNormalVector[:,1],FaceNormalVector[:,2],FaceNormalVector[:,3],lengthscale =:0.5,arrowsize=:0)