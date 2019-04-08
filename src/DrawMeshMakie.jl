function DrawMeshMakie(Val,ValName,P1,P2,P3,ColorMap)
#Draws a triangulated mesh with Makie/AbstractPlotting (08/4/19) with a set value over the meshes faces
#If you dont have a value set Val to:  Val=ones(size(P1[:,1]))
#And ValName to []; 


	#Val - value you want to Draw
	#ValName - Name of the value
	#P1 P2 P3 - Tri points in a list. Each row (of all three) representing a tri
	#ColorMap - The colormap you want to use

	#using Makie
	#using AbstractPlotting

	#Finding respective places in the to map the colours too
	MaxVal=maximum(Val);
	MinVal=minimum(Val);
	cmapsz=length(ColorMap[:,1]);
	#place around 0-1
	Val=Val.-MinVal;
	Val=Val./maximum(Val);
	Val=(Val.*(cmapsz-1)).+1
	Val=round.(Val) 
	Val=convert(Array{Int64,1},Val)

	#Applying a set color to each tri
	x=vec([P1[:,1] P2[:,1] P3[:,1]])
	y=vec([P1[:,2] P2[:,2] P3[:,2]])
	z=vec([P1[:,3] P2[:,3] P3[:,3]])
	#list in seq order (each 3 cols in list (row vec) is a tri)
	n=length(P1[:,1]);
	indices=[1:n (1:n).+n (1:n).+2*n];
	indices=indices';
	indices=vec(indices);

	#Assign into struc #RGBA color
	#Red across all tris
	color = fill(RGBAf0(1,0,0,1),length(x))
	#Green for tri no 1
	for i=1:n
		color[i]	=RGBAf0(ColorMap[Val[i],1],ColorMap[Val[i],2],ColorMap[Val[i],3],1)
		color[i+n]	=RGBAf0(ColorMap[Val[i],1],ColorMap[Val[i],2],ColorMap[Val[i],3],1)
		color[i+2*n]=RGBAf0(ColorMap[Val[i],1],ColorMap[Val[i],2],ColorMap[Val[i],3],1)
	end
	##set equal axis lims
	AxMin=minimum(vec([x y z]))
	AxMax=maximum(vec([x y z]))
	Seperation=AxMax-AxMin;
	limits = FRect3D([AxMin AxMin AxMin], [Seperation Seperation Seperation])
	p1 = Scene()
	#Draw tris (colored)
	mesh!(p1,x, y, z, indices,color=color,shading=false,limits=limits) #,limits=limits
	#Draw wireframe of mesh over the top
	wireframe!(p1[end][1], color = (:black, 0.6), linewidth = 3,limits=limits)
	#lims=p1.limits #->to check the limits are set correctly

	#Add colorbar
	colorbar = fill(RGBf0(1,0,0),length(ColorMap[:,1]))
	for i=1:cmapsz
		colorbar[i]	   =RGBf0(ColorMap[i,1],ColorMap[i,2],ColorMap[i,3])
	end
	p2 = Scene()
	cam2d!(p2)
	c = to_colormap(colorbar)
	image!(p2,
			range(0, stop = 1, length = length(c)),
			range(MinVal, stop = MaxVal, length = length(c)),
			reshape(c, (1, length(c))),
			show_axis = true)
	#Set some properties on the colourbar axis
	axis = p2[Axis]
	axis[:names, :axisnames] = ("", ValName)
	axis[:ticks, :textsize] = (0, 5)
	axis[:grid, :linewidth] = (0, 1)

	#Group the colourbar and mesh in one plot
	scene = vbox(p1, p2, sizes = [0.7, 0.3])

end