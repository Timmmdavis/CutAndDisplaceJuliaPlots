## CutAndDisplaceJuliaPlots Example<br>
<br>
using Revise<br>
using CutAndDisplaceJulia<br>
using CutAndDisplaceJuliaPlots<br>
<br>
# Correctly extract module path and run test script in CutAndDisplaceJulia to get a circular mesh<br>
ModuleDir = dirname(dirname(pathof(CutAndDisplaceJulia)))<br>
TestScript = joinpath(ModuleDir, "test", "test_TDvsEshelbyPennyCrack.jl")<br>
include(TestScript) # Run this to create a mesh we can visualize<br>
<br>
# Draw mesh with filled faces based on how much they shear<br>
CutAndDisplaceJuliaPlots.DrawMeshMakieFilledFaces(TotalShearing, "TotalShearing", P1, P2, P3, cmap2)<br>
<br>
# Draw basic mesh without filled faces<br>
fig = CutAndDisplaceJuliaPlots.DrawMeshMakie(P1, P2, P3)<br>
<br>
# Create plot with mesh boundary<br>
fig = plot(layout = 1, legend = false, xlabel = "X", ylabel = "Y", title = "Mesh Boundary")<br>
CutAndDisplaceJuliaPlots.PlotMeshBoundary(MidPoint, P1, P2, P3, FaceNormalVector, fig)<br>
