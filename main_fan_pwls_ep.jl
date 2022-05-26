using MIRT: ImageGeom, SinoGeom
using MAT
include("Reg1.jl")
include("rlalm_2d.jl")

down = 1
#sg = SinoGeom.sino_geom_ge1( ; units = :mm, strip_width = "d", down = down)
ig = image_geom( ; nx = 420, dx = 500/512)
#ig.mask = ig.circ() .> 0 immutable struct -_-

#load directories with mat files
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/kappa.mat")
kappa = vars["kappa"]
# vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/denom.mat")
# denom = vars["denom"]
# vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/sino_fan.mat")
# sino_fan = vars["sino"]
# vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/wi.mat")
# wi = vars["wi"]
# vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/xfbp.mat")
# xfbp = vars["xfbp"]
#setup edge preserving regularizer
nIter = 50
nblock = 24
l2b = 16

delta = 10
pot_arg = ["lange3"]
pot_arg2 = [delta]
R = Reg1(sqrt.(kappa), beta = 2^l2b, pot_arg = pot_arg, pot_arg2 = pot_arg2)

vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/xvalue.mat")
x = vars["x"]
