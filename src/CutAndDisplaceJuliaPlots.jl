module CutAndDisplaceJuliaPlots
# Import the packages we need
import Makie  # Use import instead of using to avoid namespace conflicts
import GLMakie  # Use import instead of using to avoid namespace conflicts
# Optionally, if you want to provide backward compatibility:
const AbstractPlotting = Makie
using DelimitedFiles
using Plots
using CutAndDisplaceJulia

# Define your own color types right in the module
struct RGBAf0
    r::Float32
    g::Float32
    b::Float32
    a::Float32
end

struct RGBf0
    r::Float32
    g::Float32
    b::Float32
end

# Constructor methods
RGBAf0(r::Number, g::Number, b::Number, a::Number) = RGBAf0(Float32(r), Float32(g), Float32(b), Float32(a))
RGBf0(r::Number, g::Number, b::Number) = RGBf0(Float32(r), Float32(g), Float32(b))

# Export these types
export RGBAf0, RGBf0

# Export these for external use
include("DrawMeshMakie.jl")
include("DrawMeshMakieFilledFaces.jl")
include("CptReader.jl")
include("LoadData2.jl")
include("SetLookDirection.jl")
include("PlotMeshBoundary.jl")

# Also export the functions
export DrawMeshMakie, DrawMeshMakieFilledFaces, CptReader, LoadData2, SetLookDirection, PlotMeshBoundary

# Load the colormap
Fid = CutAndDisplaceJuliaPlots.LoadData2(CutAndDisplaceJuliaPlots,"Ccool-warm2.cpt")
(cmap2) = CutAndDisplaceJuliaPlots.CptReader(Fid);
cmap2 = cmap2./255; # RGB2Colour
export cmap2

end # module