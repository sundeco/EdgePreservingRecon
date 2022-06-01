using MIRT: ImageGeom, SinoGeom
using MAT
using MATLAB
using Images
include("Reg1.jl")
include("rlalm_2d.jl")

down = 1
usemat = true #test by using matlab


#sg = SinoGeom.sino_geom_ge1( ; units = :mm, strip_width = "d", down = down)
ig = image_geom( ; nx = 420, dx = 500/512)
#ig.mask = ig.circ() .> 0 immutable struct -_-

#load directories with mat files
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/kappa.mat")
kappa = vars["kappa"]
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/denom.mat")
denom = vars["denom"]
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/sino_fan.mat")
sino = vars["sino"]
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/wi.mat")
wi = vars["wi"]
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/xfbp.mat")
xfbp = vars["xfbp"]
#setup edge preserving regularizer
nIter = 10
nblock = 24
l2b = 16

delta = 10
pot_arg = ["lange3"]
pot_arg2 = [delta]
R = Reg1(sqrt.(kappa), beta = 2^l2b, pot_arg = pot_arg, pot_arg2 = pot_arg2)

vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/xvalue.mat")
x = vars["x"]

vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/2d/slice420.mat")
xtrue = vars["xtrue"]
A = [1] #change this when 3d backprojection has been implemented
mask = ig.circ() .> 0
x = pwls_ep_os_rlalm_2d(xfbp[mask], A, sino, R ;
niter = nIter, wi = wi, pixmax = Inf, denom = denom, chat = false,
usemat = usemat, isave = "last", nblock = nblock, xtrue = xtrue, mask = mask)

y = MIRT.embed(x, R.mask)[:,:,1]
jim(y', clim = (800, 1200))

#assess_ssim(y, xtrue) this ssim function is wrong!
println(sqrt(sum((y-xtrue).^2/sum(R.mask)))) #RMSE gives right result
