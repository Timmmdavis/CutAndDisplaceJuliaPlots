# CutAndDisplaceJuliaPlots


Example:  
using Revise
using CutAndDisplaceJulia
using CutAndDisplaceJuliaPlots
# Correctly extract module path and run test script
ModuleDir = dirname(dirname(pathof(CutAndDisplaceJulia)))
TestScript = joinpath(ModuleDir, "test", "test_TDvsEshelbyPennyCrack.jl")
include(TestScript) #Run this to create a mesh we can visulise
CutAndDisplaceJuliaPlots.DrawMeshMakieFilledFaces(TotalShearing,"TotalShearing",P1,P2,P3,cmap2)
fig=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
fig = plot(layout = 1, legend = false, xlabel = "X", ylabel = "Y", title = "Mesh Boundary")
CutAndDisplaceJuliaPlots.PlotMeshBoundary(MidPoint,P1,P2,P3,FaceNormalVector,fig)

