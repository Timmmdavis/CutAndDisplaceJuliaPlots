function DrawMeshMakie(P1,P2,P3)
#Draws a triangulated mesh with Makie/AbstractPlotting (08/4/19) with a set value over the meshes faces
#If you dont have a value set Val to:  Val=ones(size(P1[:,1]))
#And ValName to []; 


	#Val - value you want to Draw
	#ValName - Name of the value
	#P1 P2 P3 - Tri points in a list. Each row (of all three) representing a tri
	#ColorMap - The colormap you want to use

	#using Makie
	#using AbstractPlotting

	#Applying a set color to each tri
	x=vec([P1[:,1] P2[:,1] P3[:,1]])
	y=vec([P1[:,2] P2[:,2] P3[:,2]])
	z=vec([P1[:,3] P2[:,3] P3[:,3]])
	#list in seq order (each 3 cols in list (row vec) is a tri)
	n=length(P1[:,1]);
	indices=[1:n (1:n).+n (1:n).+2*n];
	indices=indices';
	indices=vec(indices);

	##set equal axis lims
	AxMin=minimum(vec([x y z]))
	AxMax=maximum(vec([x y z]))
	Seperation=AxMax-AxMin;
	limits = FRect3D([AxMin AxMin AxMin], [Seperation Seperation Seperation])
	p1 = Scene()
	#Draw tris (colored)
	mesh!(p1,x, y, z, indices,shading=false,limits=limits, color = :blue) #,limits=limits
	#Draw wireframe of mesh over the top
	wireframe!(p1[end][1], color = (:black, 0.6), linewidth = 3,limits=limits)
	#lims=p1.limits #->to check the limits are set correctly

	return p1

end