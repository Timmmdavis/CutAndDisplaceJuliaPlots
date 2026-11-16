## CutAndDisplaceJuliaPlots Example<br>
<br>
using Revise<br>
using Plots<br>
using UnicodePlots<br>
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
<br>
## Crack views (fast, for propagation runs)<br>
<br>
Square figure with three aligned views of the opened crack in its own frame (face-on, edge-on along strike, edge-on down-dip), triangles coloured by opening and the tip line coloured by KI/Kc (red >1 grows, blue <1). Drawn with CairoMakie (no GPU), about 0.2 s per frame.<br>
<br>
fig = CutAndDisplaceJuliaPlots.PlotCrackState(P1, P2, P3, Dn; TipXYZ = [XMid YMid ZMid], TipK = StrainEnergyV, Kc = KCrit, meta = (step = i,))<br>
CutAndDisplaceJuliaPlots.SaveFig("step.png", fig)<br>
<br>
# Every 5th step of a finished run (reads its *-MeshWithDisp-ForDrawing-*.txt and *-MeshEdges-ForDrawing-*.txt files)<br>
CutAndDisplaceJuliaPlots.RenderRun(raw"path\to\run"; every = 5)<br>
<br>
# Check the layout on real, rotated and twisted cracks (writes path\to\run\FigureTest)<br>
CutAndDisplaceJuliaPlots.TestFigure(raw"path\to\run")<br>
