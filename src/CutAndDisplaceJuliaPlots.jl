module CutAndDisplaceJuliaPlots
using Makie
using AbstractPlotting
using DelimitedFiles

include("DrawMeshMakie.jl")
include("DrawMeshMakieFilledFaces.jl")
include("CptReader.jl")
include("LoadData2.jl")
include("SetLookDirection.jl")

#Load the colormap
Fid=CutAndDisplaceJuliaPlots.LoadData2(CutAndDisplaceJuliaPlots,"Ccool-warm2.cpt")
(cmap2)=CutAndDisplaceJuliaPlots.CptReader(Fid);
cmap2=cmap2./255; #RGB2Colour
export cmap2

end # module


