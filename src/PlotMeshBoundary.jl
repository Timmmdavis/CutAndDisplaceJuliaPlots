function PlotMeshBoundary(MidPoint,P1,P2,P3,FaceNormalVector,fig)



(UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint)=
CutAndDisplaceJulia.GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)

(SortedTriangles,ConnectedEdge)=CutAndDisplaceJulia.ConnectedConstraints(P1,P2,P3,MidPoint);



LeadingPoint=[0. 0. 0.]
TrailingPoint=[0. 0. 0.]
InnerPoint   =[0. 0. 0.]    
LeadingPointOld  =[NaN NaN NaN] 
TrailingPointOld =[NaN NaN NaN] 
InnerPointOld    =[NaN NaN NaN] 
BackPoint    =[NaN NaN NaN] 


idxa=1
idxb=3

lps=0
#For each edge loop
for i=1:length(UniqueEdges)

  b=vec(UniqueEdges[i])

  for j=b

    

    #Extract the points on the current bit of the edge
    (LeadingPoint,~) =GrabPointNew2(LeadingPoints,P1,P2,P3,j)
    (TrailingPoint,~)=GrabPointNew2(TrailingPoints,P1,P2,P3,j)
    (InnerPoint,Idx) =GrabPointNew2(InnerPoints,P1,P2,P3,j)

    #lps+=1
    #if lps==1
    #  plot([LeadingPoint],[TrailingPoint])
    #else
      Plots.plot!([LeadingPoint[idxa],TrailingPoint[idxa]],[LeadingPoint[idxb],TrailingPoint[idxb]],c=(:red), lab="")
    #end
    
    end

end 

return fig

end


function GrabPointNew2(PointsIdxList,P1,P2,P3,j)
#Extract the points on the current bit of the edge
Point=[0. 0. 0.]
Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        Point=P1[Idx,:]
        Indx=Idx
    elseif k==2
        Point=P2[Idx,:]
        Indx=Idx
    elseif k==3
        Point=P3[Idx,:]
        Indx=Idx
    end
end

return Point,Indx
end

