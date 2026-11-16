# Fast per-step crack figures (Makie, drawn with CairoMakie - no GPU needed).
#
#   fig = PlotCrackState(P1,P2,P3,Dn; TipXYZ, TipK, Kc, meta=(step=i, V=..., P=..., ρ=...))
#   SaveFig("step.png", fig)
#   RenderRun(dir; every=5)     # frames for a whole run from its *-MeshWithDisp-ForDrawing-*.txt files
#   TestFigure(dir)             # real, rotated and twisted cracks -> dir/FigureTest
#
# Square figure, one metres-per-pixel scale, world axes (third-angle layout):
#   x–z (looking north)                  | y–z (looking west), shares z
#   x–y (map view from above), shares x  | colour bars: opening Dn, tip KI/Kc
# The opened crack is drawn as in CutAndDisplaceJulia.ExportCrackMesh ("MeshFilled"): walls at
# ±Exag*Dn/2, triangles coloured by Dn; opening exaggerated (across the crack only) so the widest
# aperture is at least 1/10 of the plot. Tip line coloured by KI/Kc on all views (red >1, blue <1).

_figure(w, h) = try MK.Figure(size = (w, h)) catch; MK.Figure(resolution = (w, h)) end
_noshade() = isdefined(MK, :NoShading) ? MK.NoShading : false
SaveFig(path, fig) = try MK.save(path, fig; px_per_unit = 1) catch; MK.save(path, fig) end

# Best-fit plane of the crack: n̂ normal, û along strike (horizontal), v̂ up-dip
function CrackFrame(P1, P2, P3)
    V = vcat(P1, P2, P3)
    c = vec(mean(V, dims = 1))
    n̂ = svd(V .- c').V[:, 3]
    ẑ = [0.0, 0.0, 1.0]
    û = abs(dot(n̂, ẑ)) > 0.99 ? [1.0, 0.0, 0.0] : normalize(cross(ẑ, n̂))
    v̂ = cross(n̂, û)
    v̂[3] < 0 && (v̂ = -v̂)
    return c, û, v̂, n̂
end

# Opening at mesh vertices: mean of the adjacent faces, zero on the crack edge (tips close)
function VertexDn(P1, P2, P3, Dn)
    n = size(P1, 1)
    key = Dict{NTuple{3,Float64},Int}()
    ids = zeros(Int, n, 3)
    for (k, P) in enumerate((P1, P2, P3)), i in 1:n
        t = (round(P[i, 1], digits = 6), round(P[i, 2], digits = 6), round(P[i, 3], digits = 6))
        ids[i, k] = get!(key, t, length(key) + 1)
    end
    nv = length(key)
    sm = zeros(nv); cnt = zeros(nv)
    for k in 1:3, i in 1:n
        sm[ids[i, k]] += Dn[i]; cnt[ids[i, k]] += 1
    end
    vd = sm ./ cnt
    ec = Dict{Tuple{Int,Int},Int}()
    for i in 1:n, (p, q) in ((1, 2), (2, 3), (3, 1))
        e = minmax(ids[i, p], ids[i, q]); ec[e] = get(ec, e, 0) + 1
    end
    for (e, m) in ec
        m == 1 && (vd[e[1]] = 0.0; vd[e[2]] = 0.0)
    end
    return ids, vd
end

# Crack walls as in CutAndDisplaceJulia.ExportCrackMesh ("MeshFilled"): every vertex moved by
# ±E*Dn/2 along its averaged normal, Dn averaged from adjacent faces and zero on the crack edge.
function CrackWalls(P1, P2, P3, Dn, E, n̂)
    n = size(P1, 1)
    ids, vd = VertexDn(P1, P2, P3, Dn)
    nv = maximum(ids)
    X = zeros(nv, 3)
    for (k, P) in enumerate((P1, P2, P3)), i in 1:n
        X[ids[i, k], :] .= P[i, :]
    end
    Nv = zeros(nv, 3)
    for i in 1:n
        fn = cross(P2[i, :] .- P1[i, :], P3[i, :] .- P1[i, :])
        nf = norm(fn); nf == 0 && continue
        fn ./= nf; dot(fn, n̂) < 0 && (fn .*= -1)
        for k in 1:3
            Nv[ids[i, k], :] .+= fn
        end
    end
    for j in 1:nv
        nn = norm(Nv[j, :]); Nv[j, :] .= nn > 0 ? Nv[j, :] ./ nn : n̂
    end
    W1 = X .+ (E .* vd ./ 2) .* Nv
    W2 = X .- (E .* vd ./ 2) .* Nv
    return ids, X, W1, W2
end

# Crack-aligned viewing frame: û along strike (azimuth in [0,180)), v̂ up-dip, n̂ normal (right-handed)
function ViewFrame(P1, P2, P3)
    c, û, v̂, n̂ = CrackFrame(P1, P2, P3)
    az = mod(atand(û[1], û[2]), 360)
    if az >= 180
        û = -û; n̂ = -n̂; az -= 180
    end
    v̂ = cross(n̂, û)
    dip = acosd(clamp(abs(n̂[3]), 0, 1))
    n0 = dot(c, n̂)                                            # plane offset: across-crack axis reads ± about 0
    return (û = û, v̂ = v̂, n̂ = n̂, az = az, dip = dip, n0 = n0)
end

# Boundary (tip) edges of the mesh as vertex pairs
function BoundaryEdges(ids)
    ec = Dict{Tuple{Int,Int},Int}()
    for i in 1:size(ids, 1), (p, q) in ((1, 2), (2, 3), (3, 1))
        e = minmax(ids[i, p], ids[i, q]); ec[e] = get(ec, e, 0) + 1
    end
    return [e for (e, m) in ec if m == 1]
end

# One orthographic view of both walls. Pr1/Pr2: wall vertices in view coordinates (nv×3);
# plot columns (ia, ib), depth column id, sgn=+1 viewer at +∞ of id.
function _drawview!(ax, ids, Pr1, Pr2, Dn, ia, ib, id, sgn, Dmax; edgealpha = 0.18)
    n = size(ids, 1)
    tris = [(W, i) for W in (Pr1, Pr2) for i in 1:n]
    depth = [sgn * mean(W[ids[i, :], id]) for (W, i) in tris]
    order = sortperm(depth)                                   # far first (painter's algorithm)
    verts = zeros(3 * length(order), 2); cols = zeros(3 * length(order))
    ex = Float64[]; ey = Float64[]
    for (m, o) in enumerate(order)
        (W, i) = tris[o]
        for k in 1:3
            verts[3(m - 1) + k, 1] = W[ids[i, k], ia]
            verts[3(m - 1) + k, 2] = W[ids[i, k], ib]
            cols[3(m - 1) + k] = Dn[i]
        end
        p = W[ids[i, :], :]
        append!(ex, [p[1, ia], p[2, ia], p[2, ia], p[3, ia], p[3, ia], p[1, ia]])
        append!(ey, [p[1, ib], p[2, ib], p[2, ib], p[3, ib], p[3, ib], p[1, ib]])
    end
    faces = collect(reshape(collect(1:size(verts, 1)), 3, :)')
    mp = MK.mesh!(ax, verts, faces; color = cols, colormap = :viridis, colorrange = (0, Dmax),
                  shading = _noshade())
    MK.linesegments!(ax, ex, ey, color = (:black, edgealpha), linewidth = 0.4)
    return mp
end

"""
PlotCrackState(P1,P2,P3,Dn; TipXYZ, TipK, Kc, meta, Dnmax, Exag, lims, pxmax)

Three aligned orthographic views of the opened crack in world coordinates:
  x–z (looking north, +y)            | y–z (looking west, -x), shares z
  x–y (map view from above), shares x | colour bars: opening Dn, tip KI/Kc
Walls at ±Exag*Dn/2 along the crack normal (as ExportCrackMesh), triangles coloured by opening,
tip line coloured by KI/Kc on all views. One metres-per-pixel scale for every panel.
lims = (xlo,xhi,ylo,yhi,zlo,zhi) fixes the views (RenderRun sets it for a whole run).
"""
function PlotCrackState(P1, P2, P3, Dn; TipXYZ = zeros(0, 3), TipK = Float64[], Kc = 5e7,
                        meta = NamedTuple(), Dnmax = nothing, Exag = nothing, lims = nothing, pxmax = 800)
    n = size(P1, 1)
    F = ViewFrame(P1, P2, P3)                                  # only for the crack normal and strike/dip text
    V = vcat(P1, P2, P3)
    Dmax = Dnmax === nothing ? max(maximum(Dn), eps()) : Dnmax
    if Exag === nothing
        span = max(maximum(V[:, 1]) - minimum(V[:, 1]), maximum(V[:, 2]) - minimum(V[:, 2]),
                   maximum(V[:, 3]) - minimum(V[:, 3]))
        Exag = max(1.0, 0.1 * 1.1 * span / Dmax)              # widest aperture ≥ 1/10 of the plot
    end
    ids, X, W1, W2 = CrackWalls(P1, P2, P3, Dn, Exag, F.n̂)
    if lims === nothing
        A = vcat(W1, W2)
        lo = vec(minimum(A, dims = 1)); hi = vec(maximum(A, dims = 1))
        pad = 0.05 * maximum(hi .- lo)
        lims = (lo[1] - pad, hi[1] + pad, lo[2] - pad, hi[2] + pad, lo[3] - pad, hi[3] + pad)
    end
    (xlo, xhi, ylo, yhi, zlo, zhi) = lims
    # Square drawing: x–z box made square (side L); y range at least 0.32 L (room for colour bars)
    L = max(xhi - xlo, zhi - zlo)
    xc = (xlo + xhi) / 2; zc = (zlo + zhi) / 2
    xlo, xhi = xc - L / 2, xc + L / 2
    zlo, zhi = zc - L / 2, zc + L / 2
    yr = max(yhi - ylo, 0.32 * L)
    yc = (ylo + yhi) / 2; ylo, yhi = yc - yr / 2, yc + yr / 2
    s = pxmax / (L + yr)                                      # pixels per metre, same in every panel

    # ---- text
    tipz = maximum(V[:, 3])
    parts = String[]
    push!(parts, @sprintf("top tip z = %.1f m", tipz))
    haskey(meta, :V) && isfinite(meta.V) && push!(parts, @sprintf("fluid volume = %.4g m³", meta.V))
    haskey(meta, :P) && isfinite(meta.P) && push!(parts, @sprintf("P(ref) = %.3g MPa", meta.P / 1e6))
    haskey(meta, :ρ) && isfinite(meta.ρ) && push!(parts, @sprintf("ρ̄ fluid = %.1f kg/m³", meta.ρ))
    push!(parts, @sprintf("max Dn = %.3g m", maximum(Dn)))
    push!(parts, "$n elements")
    ttl = (haskey(meta, :step) ? "Step $(meta.step)" : "Crack") * "  |  " * join(parts, "  |  ")
    etxt = Exag > 1.5 ? @sprintf("opening ×%.0f", Exag) : "opening ×1 (true scale)"
    sub = @sprintf("Best-fit crack plane: strike %.0f°, dip %.0f°.  Opening exaggerated ×%.0f along the crack normal only.",
                   F.az, F.dip, Exag)

    # ---- layout: exact pixel sizes per axis, figure resized to fit
    fig = _figure(1400, 1400)
    MK.Label(fig[1, 1:2], ttl, fontsize = 16, tellwidth = false)
    MK.Label(fig[2, 1:2], sub, fontsize = 13, color = :grey30, tellwidth = false)
    axF = MK.Axis(fig[3, 1], width = s * L, height = s * L, xlabel = "x [m]", ylabel = "z [m]",
                  title = "x–z (looking north, +y)")
    axS = MK.Axis(fig[3, 2], width = s * yr, height = s * L, xlabel = "y [m]", ylabel = "z [m]",
                  title = "y–z (looking west, −x)")
    axM = MK.Axis(fig[4, 1], width = s * L, height = s * yr, xlabel = "x [m]", ylabel = "y [m]",
                  title = "x–y (map view, from above)")
    mp = _drawview!(axF, ids, W1, W2, Dn, 1, 3, 2, -1, Dmax)  # viewer at -y
    _drawview!(axS, ids, W1, W2, Dn, 2, 3, 1, 1, Dmax; edgealpha = 0.1)   # viewer at +x
    _drawview!(axM, ids, W1, W2, Dn, 1, 2, 3, 1, Dmax; edgealpha = 0.1)   # viewer above
    if zhi >= 0
        for ax in (axF, axS)
            MK.hlines!(ax, [0.0], color = :saddlebrown, linestyle = :dash, linewidth = 2)
        end
    end

    # ---- tip line coloured by K/Kc, all views
    BE = BoundaryEdges(ids)
    kseg = fill(NaN, length(BE))
    good = size(TipXYZ, 1) > 0 ? findall(i -> all(isfinite, TipXYZ[i, :]) && isfinite(TipK[i]), 1:size(TipXYZ, 1)) : Int[]
    if !isempty(good)
        T3 = TipXYZ[good, :]; Kg = TipK[good] ./ Kc
        for (q, (a, b)) in enumerate(BE)
            mid = (X[a, :] .+ X[b, :]) ./ 2; len = norm(X[a, :] .- X[b, :])
            d2 = vec(sum((T3 .- mid') .^ 2, dims = 2))
            j = argmin(d2)
            sqrt(d2[j]) < 0.75 * len && (kseg[q] = Kg[j])
        end
    end
    kcmap = MK.cgrad([:navy, :deepskyblue, :orangered, :darkred], [0.0, 0.499, 0.501, 1.0])
    function tipline!(ax, ia, ib; lw = 7)
        tx = Float64[]; ty = Float64[]; tc = Float64[]; gx = Float64[]; gy = Float64[]
        for (q, (a, b)) in enumerate(BE)
            if isfinite(kseg[q])
                append!(tx, [X[a, ia], X[b, ia]]); append!(ty, [X[a, ib], X[b, ib]]); append!(tc, [kseg[q], kseg[q]])
            else
                append!(gx, [X[a, ia], X[b, ia]]); append!(gy, [X[a, ib], X[b, ib]])
            end
        end
        isempty(tx) && isempty(gx) || MK.linesegments!(ax, vcat(tx, gx), vcat(ty, gy), color = :black, linewidth = lw + 3)
        isempty(gx) || MK.linesegments!(ax, gx, gy, color = :grey70, linewidth = lw)
        return isempty(tx) ? nothing :
               MK.linesegments!(ax, tx, ty; color = clamp.(tc, 0, 2), colormap = kcmap, colorrange = (0, 2), linewidth = lw)
    end
    kp = tipline!(axF, 1, 3)
    tipline!(axS, 2, 3; lw = 5)
    tipline!(axM, 1, 2; lw = 5)

    MK.limits!(axF, xlo, xhi, zlo, zhi)
    MK.limits!(axS, ylo, yhi, zlo, zhi)
    MK.limits!(axM, xlo, xhi, ylo, yhi)
    MK.linkyaxes!(axF, axS); MK.linkxaxes!(axF, axM)

    # exaggeration marker: double-headed arrow along the crack normal as seen in each view
    # (drawn only where the normal is not end-on), in the panel's lower-left corner
    for (ax, ia, ib, alo, ahi, blo, bhi) in ((axF, 1, 3, xlo, xhi, zlo, zhi), (axS, 2, 3, ylo, yhi, zlo, zhi),
                                            (axM, 1, 2, xlo, xhi, ylo, yhi))
        d = [F.n̂[ia], F.n̂[ib]]
        norm(d) < 0.3 && continue
        d ./= norm(d)
        half = 0.1 * min(ahi - alo, bhi - blo)
        cx = alo + 0.06 * (ahi - alo) + half * abs(d[1]); cy = blo + 0.06 * (bhi - blo) + half * abs(d[2])
        p0 = (cx - half * d[1], cy - half * d[2]); p1 = (cx + half * d[1], cy + half * d[2])
        MK.lines!(ax, [p0[1], p1[1]], [p0[2], p1[2]], color = :black, linewidth = 1.5)
        MK.scatter!(ax, [p0[1], p1[1]], [p0[2], p1[2]], color = :black, markersize = 9, marker = :circle)
        MK.text!(ax, cx, cy, text = etxt, align = (:left, :bottom), fontsize = 11, offset = (6, 6))
    end

    # colour bars, horizontal, in the empty bottom-right corner
    cb = MK.GridLayout(fig[4, 2], tellwidth = false, tellheight = false)
    MK.Colorbar(cb[1, 1], mp, label = "Opening Dn [m]", vertical = false, flipaxis = false,
                width = MK.Relative(0.9), height = 12, labelsize = 12, ticklabelsize = 10)
    kp === nothing || MK.Colorbar(cb[2, 1], kp, label = "Tip line KI/Kc: red >1 grows, blue <1 (grey: no value)",
                                  vertical = false, flipaxis = false, width = MK.Relative(0.9), height = 12,
                                  labelsize = 12, ticklabelsize = 10, ticks = [0, 0.5, 1, 1.5, 2])
    try
        MK.resize_to_layout!(fig)
    catch
    end
    return fig
end

# Quick check of the figure on real steps and on rotated / twisted copies of the crack.
# Writes PNGs to <dir>/FigureTest and prints their paths.
function TestFigure(dir)
    steps = FindSteps(dir)
    out = joinpath(dir, "FigureTest"); mkpath(out)
    files = String[]
    for st in (steps[cld(length(steps), 3)], steps[end])
        d = ReadStep(st)
        f = PlotCrackState(d.P1, d.P2, d.P3, d.Dn; TipXYZ = d.TipXYZ, TipK = d.TipK, Kc = d.Kc, meta = (step = st.i,))
        push!(files, joinpath(out, "step_$(st.i).png")); SaveFig(files[end], f)
    end
    d = ReadStep(steps[cld(length(steps), 2)])
    c = vec(mean(vcat(d.P1, d.P2, d.P3), dims = 1))
    # rotated: strike 30°, dip 60°
    R = [cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1] * [1 0 0; 0 cosd(30) -sind(30); 0 sind(30) cosd(30)]
    rot(P) = ((P .- c') * R') .+ c'
    f = PlotCrackState(rot(d.P1), rot(d.P2), rot(d.P3), d.Dn; TipXYZ = rot(d.TipXYZ), TipK = d.TipK, Kc = d.Kc,
                       meta = (step = "rotated",))
    push!(files, joinpath(out, "rotated.png")); SaveFig(files[end], f)
    # twisted: rotation about the vertical axis increasing with height (40° over the crack height)
    zr = extrema(vcat(d.P1, d.P2, d.P3)[:, 3])
    function twist(P)
        Q = copy(P)
        for i in 1:size(P, 1)
            θ = 40 * (P[i, 3] - zr[1]) / (zr[2] - zr[1])
            x = P[i, 1] - c[1]; y = P[i, 2] - c[2]
            Q[i, 1] = c[1] + cosd(θ) * x - sind(θ) * y; Q[i, 2] = c[2] + sind(θ) * x + cosd(θ) * y
        end
        Q
    end
    f = PlotCrackState(twist(d.P1), twist(d.P2), twist(d.P3), d.Dn; TipXYZ = twist(d.TipXYZ), TipK = d.TipK, Kc = d.Kc,
                       meta = (step = "twisted",))
    push!(files, joinpath(out, "twisted.png")); SaveFig(files[end], f)
    println.(files)
    return files
end

# ------------------------------------------------------------------ reading a run's output files
function _readtable(f)
    rows = Vector{Vector{Float64}}()
    for (k, l) in enumerate(eachline(f))
        k == 1 && continue
        s = split(strip(l))
        isempty(s) && continue
        push!(rows, parse.(Float64, s))
    end
    return isempty(rows) ? zeros(0, 0) : reduce(vcat, permutedims.(rows))
end

"Steps of the latest run in `dir` (files `i-p-MeshWithDisp-ForDrawing-Rand.txt`)."
function FindSteps(dir)
    re = r"^(\d+)-(\d+)-MeshWithDisp-ForDrawing-(\d+)\.txt$"
    hits = Tuple{Int,String,String,String}[]
    for f in readdir(dir)
        m = match(re, f)
        m === nothing || push!(hits, (parse(Int, m[1]), m[2], m[3], f))
    end
    isempty(hits) && error("no MeshWithDisp-ForDrawing files in $dir")
    rands = [h[3] for h in hits]
    R = argmax(r -> count(==(r), rands), unique(rands))
    hits = sort(filter(h -> h[3] == R, hits), by = first)
    return [(i = h[1], mesh = joinpath(dir, h[4]),
             edges = joinpath(dir, "$(h[1])-$(h[2])-MeshEdges-ForDrawing-$(h[3]).txt")) for h in hits]
end

function ReadStep(st)
    M = _readtable(st.mesh)
    P1 = M[:, 1:3]; P2 = M[:, 4:6]; P3 = M[:, 7:9]; Dn = M[:, 10]
    E = isfile(st.edges) ? _readtable(st.edges) : zeros(0, 5)
    TipXYZ = size(E, 1) > 0 ? E[:, 1:3] : zeros(0, 3)
    TipK = size(E, 1) > 0 ? E[:, 4] : Float64[]
    Kc = size(E, 1) > 0 ? E[1, 5] : 5e7
    return (P1 = P1, P2 = P2, P3 = P3, Dn = Dn, TipXYZ = TipXYZ, TipK = TipK, Kc = Kc)
end

function _readlog(dir)
    f = joinpath(dir, "SolverLog.csv")
    isfile(f) || return Dict{Int,NamedTuple}()
    (d, h) = readdlm(f, ',', header = true)
    h = vec(String.(h))
    col(name) = findfirst(==(name), h)
    out = Dict{Int,NamedTuple}()
    for r in eachrow(d)
        it = Int(r[col("Iter")])
        haskey(out, it) && continue                   # first row of a step is the main solver
        P = r[col("P")]; s = r[col("s")]
        out[it] = (V = Float64(r[col("OpenedVolume")]), P = P isa Number ? Float64(P) : NaN,
                   s = s isa Number ? Float64(s) : NaN)
    end
    return out
end

"""
RenderRun(dir; every=1, outdir=joinpath(dir,"frames"), ρ0=nothing)

Renders every `every`-th step of the latest run in `dir` with fixed colour/axis limits.
Pass ρ0 (fluid reference density) to show mean density = s*ρ0 for gas runs.
"""
function RenderRun(dir; every = 1, outdir = joinpath(dir, "frames"), ρ0 = nothing)
    steps = FindSteps(dir)[1:every:end]
    data = [ReadStep(st) for st in steps]
    log = _readlog(dir)
    Dmax = maximum(maximum(d.Dn) for d in data)
    Q = vcat([vcat(d.P1, d.P2, d.P3) for d in data]...)
    lo = vec(minimum(Q, dims = 1)); hi = vec(maximum(Q, dims = 1))
    hi[3] = max(hi[3], 0.0)                                  # keep the surface in view
    span = maximum(begin W = vcat(d.P1, d.P2, d.P3); maximum(maximum(W, dims = 1) .- minimum(W, dims = 1)) end for d in data)
    E = max(1.0, 0.1 * 1.1 * span / Dmax)                    # same exaggeration in every frame
    pad = 0.04 * max(hi[1] - lo[1], hi[3] - lo[3])
    hw = 0.6 * E * Dmax
    L = (lo[1] - pad, hi[1] + pad, lo[2] - hw - pad, hi[2] + hw + pad, lo[3] - pad, hi[3] + pad)
    mkpath(outdir)
    times = Float64[]
    for (st, d) in zip(steps, data)
        lg = get(log, st.i, (V = NaN, P = NaN, s = NaN))
        meta = (step = st.i, V = lg.V, P = lg.P, ρ = ρ0 === nothing ? NaN : lg.s * ρ0)
        t = @elapsed begin
            fig = PlotCrackState(d.P1, d.P2, d.P3, d.Dn; TipXYZ = d.TipXYZ, TipK = d.TipK, Kc = d.Kc,
                                 meta = meta, Dnmax = Dmax, Exag = E, lims = L)
            SaveFig(joinpath(outdir, @sprintf("frame_%04d.png", st.i)), fig)
        end
        push!(times, t)
    end
    println("Rendered $(length(steps)) frames to $outdir; first $(round(times[1], digits=2)) s, median $(round(median(times), digits=3)) s")
    return times
end
