using MIRT: ImageGeom, SinoGeom
using MAT
include("Reg1.jl")

down = 1
sg = sino_geom_ge1( ; units = :mm, strip_width = "d", down = down)
ig = image_geom( ; nx = 420, dx = 500/512, down = down)
