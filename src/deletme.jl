#Some example plotting (Can do during debug)
scene=CutAndDisplaceJuliaPlots.DrawMeshMakieFilledFaces(Dn,"Dn",P1,P2,P3,cmap2)

using Makie

#Some plotting of the mesh as we go. 
#Show current mesh:
CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
#Plots to current scene if a scene is not set by default
scatter!(Px, Py, Pz) #markersize = 50


#Show normals - scale length scale by max sep of tris
arrows!(MidPoint[:,1],MidPoint[:,2],MidPoint[:,3],
	FaceNormalVector[:,1],FaceNormalVector[:,2],FaceNormalVector[:,3],arrowsize=:0) #lengthscale =:150

#Show computed edge points:
scatter!(FeP1P2S.FeMd[:,1],FeP1P2S.FeMd[:,2],FeP1P2S.FeMd[:,3],color = :red) #,markersize  =:150
scatter!(FeP1P3S.FeMd[:,1],FeP1P3S.FeMd[:,2],FeP1P3S.FeMd[:,3],color = :green)
scatter!(FeP2P3S.FeMd[:,1],FeP2P3S.FeMd[:,2],FeP2P3S.FeMd[:,3],color = :blue)

#Show computed mid2edge vectors:
arrows!(FeP1P2S.FeMd[:,1],FeP1P2S.FeMd[:,2],FeP1P2S.FeMd[:,3],
	FeP1P2S.FeM2Ev[:,1],FeP1P2S.FeM2Ev[:,2],FeP1P2S.FeM2Ev[:,3], linecolor  = :red) #,lengthscale =:150
arrows!(FeP1P3S.FeMd[:,1],FeP1P3S.FeMd[:,2],FeP1P3S.FeMd[:,3],
	FeP1P3S.FeM2Ev[:,1],FeP1P3S.FeM2Ev[:,2],FeP1P3S.FeM2Ev[:,3], linecolor  = :green)
arrows!(FeP2P3S.FeMd[:,1],FeP2P3S.FeMd[:,2],FeP2P3S.FeMd[:,3],
	FeP2P3S.FeM2Ev[:,1],FeP2P3S.FeM2Ev[:,2],FeP2P3S.FeM2Ev[:,3] ,linecolor  = :blue)

#Show computed along edge vectors:
arrows!(FeP1P2S.FeMd[:,1],FeP1P2S.FeMd[:,2],FeP1P2S.FeMd[:,3],
	FeP1P2S.FeEv[:,1],FeP1P2S.FeEv[:,2],FeP1P2S.FeEv[:,3],arrowsize=:0.05, linecolor  = :red) #,lengthscale =:150
arrows!(FeP1P3S.FeMd[:,1],FeP1P3S.FeMd[:,2],FeP1P3S.FeMd[:,3],
	FeP1P3S.FeEv[:,1],FeP1P3S.FeEv[:,2],FeP1P3S.FeEv[:,3],arrowsize=:0.05, linecolor  = :green)
arrows!(FeP2P3S.FeMd[:,1],FeP2P3S.FeMd[:,2],FeP2P3S.FeMd[:,3],
	FeP2P3S.FeEv[:,1],FeP2P3S.FeEv[:,2],FeP2P3S.FeEv[:,3],arrowsize=:0.05,linecolor  = :blue,linewidth=:4)




using Plots

Px=[p1[:,1]; p2[:,1]; p3[:,1]];
Py=[p1[:,2]; p2[:,2]; p3[:,2]];
Pz=[p1[:,3]; p2[:,3]; p3[:,3]];
Ang=[Ang1; Ang2; Ang3]
a=2; #plot for x
b=3; #plot for z
XMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,a]
	  FeP1P3S.FeMd[FeP1P3S.FreeFlg,a]
	  FeP2P3S.FeMd[FeP2P3S.FreeFlg,a]]
YMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,b]
	  FeP1P3S.FeMd[FeP1P3S.FreeFlg,b]
	  FeP2P3S.FeMd[FeP2P3S.FreeFlg,b]]
XDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,a]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,a]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,a]]
YDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,b]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,b]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,b]]
XDir2=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,a]
	   FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,a]
	   FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,a]]
YDir2=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,b]
	   FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,b]
	   FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,b]]	
K2=[FeP1P2S.K2[FeP1P2S.FreeFlg,1]
    FeP1P3S.K2[FeP1P3S.FreeFlg,1]
    FeP2P3S.K2[FeP2P3S.FreeFlg,1]]	
K1=[FeP1P2S.K1[FeP1P2S.FreeFlg,1]
    FeP1P3S.K1[FeP1P3S.FreeFlg,1]
    FeP2P3S.K1[FeP2P3S.FreeFlg,1]]	    	     


#Functionize and return final fig (all 3) - or just group before?
gr()
scl=0.1;#length of vectors
fig = plot(reuse=false,legend=false)
for i=1:length(XMid)

	x=XMid[i]; 
	y=YMid[i]; 
	scatter!([x],[y],ms=K2[i]*6)
	xnew=x+XDir[i]*scl; 
	ynew=y+YDir[i]*scl; 
	xnew2=x+XDir2[i]*scl; 
	ynew2=y+YDir2[i]*scl; 
	xd1=[x xnew]
	yd1=[y ynew]
	xd2=[x xnew2]
	yd2=[y ynew2]
	plot!(vec(xd1),vec(yd1))
	plot!(vec(xd2),vec(yd2))

end 	

using Plots
gr()
scl=0.1;#length of vectors
fig2 = plot(reuse=false,legend=false)
scatter!(Px,Pz,ms=-Py*4) #Looking towards pos y