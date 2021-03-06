using MIRT: ImageGeom, SinoGeom
using MAT
using MATLAB
using Images
using FileIO, JLD2

include("Reg1.jl")
include("rlalm_3d.jl")
include("reshaper.jl")
include("/Users/jasonhu/Documents/julia_files/sf_from_scratch/gblock.jl")

down = 1
usemat = false

ig = image_geom( ; nx = 420, dx = 500/512, nz = 96, dz = 0.625)
cg = ct_geom(:ge1)

if usemat
    vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/kappa.mat")
    kappa = vars["kappa"]
    vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/denom.mat")
    denom = vars["denom"]
    vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/sino_cone.mat")
    sino = vars["sino"]
    vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/wi.mat")
    wi = vars["wi"]
else
    sino = FileIO.load("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/sino_true.jld2", "sino_true");
    kappa = FileIO.load("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/kappa.jld2", "kappa");
    denom = FileIO.load("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/denom.jld2", "denom");
    wi = FileIO.load("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/wi.jld2", "wi");
end
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/xfdk.mat")
xfdk = vars["xfdk"]

nIter = 2
nblock = 24
l2b = 14.5
delta = 10

vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/xtrue_crop17-80.mat")
xtrue = vars["xtrue"]

pot_arg = ["lange3"]
pot_arg2 = [delta]
mask = ig.circ() .> 0

b1 = ig.dx^2
b2 = 1/(ig.dx^2 + ig.dy^2)
b3 = 1/ig.dz^2
b4 = 1/(ig.dx^2 + ig.dz^2)
b5 = 1/(ig.dx^2+ig.dy^2+ig.dz^2)
beta = 2^l2b*[b1, b1, b2, b2, b5, b4, b5, b4, b3, b4, b5, b4, b5]
R = Reg1(sqrt.(kappa), beta = beta, offsets = "3d:26", pot_arg = pot_arg,
pot_arg2 = pot_arg2, distance_power = 0, mask = mask)

vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/3d/xvalue.mat")
x = vars["x"]
#error("a")
#denom = denom[R.mask]
A = gblock(cg, ig, nblock)
x = pwls_ep_os_rlalm_3d(xfdk[mask], A, reshaper(sino, "2d"), R ;
wi = reshaper(wi, "2d"), niter = nIter, denom = denom, chat = false,
xtrue = xtrue, mask = mask, usemat = usemat, nblock = nblock)

y = MIRT.embed(x, R.mask)[:,:,:,1]
#jim(y', clim = (800, 1200))
jim(permutedims(y, (2,1,3)), clim = (800, 1200))
