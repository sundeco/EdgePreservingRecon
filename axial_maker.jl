using MATLAB
using MIRT: ImageGeom, SinoGeom
using FileIO, JLD2
include("/Users/jasonhu/Documents/julia_files/sf_from_scratch/gcone.jl")

IOs = 10000

cg = ct_geom(:ge1)
vars = matread("/Users/jasonhu/Documents/GitHub/EdgePreservingRecon/data/phantom_crop154.mat")
phantom = vars["phantom"]
phantom = phantom[:,:,1:96]

mm2HU = 1000 / 0.02
ig_hi = image_geom( ; nx = 840, dx = 500/1024, nz = 96, dz = 0.625)
A_hi = gcone_map(cg, ig_hi)
# sino_true = A_hi * phantom #this will probably take a very long time to run
# FileIO.save("sino_true.jld2", "sino_true", sino_true)
FileIO.load("sino_true.jld2", "sino_true");

#add noise to get wi and yi
yi = poisson(IOs * exp.(-sino_true / mm2HU), factor=0.4)
var = 5
ye = var * randn(size(yi))
k = 1
zi = k * yi + ye
myerror = 1/10
zi = max.(zi, myerror)
sino = -log.(zi /(k*IOs)) * mm2HU
wi = (zi.^2)./(k*zi .+ var^2)
FileIO.save("wi.jld2", "wi", wi)

#set up target geometry
ig = image_geom( ; nx = 420, dx = 500/512, nz = 96, dz = 0.625)
A = gcone_map(cg, ig)

#set up kappa
kappa = sqrt.(div0.(A' * wi, A' * ones(size(wi))))
kappa = max.(kappa, 0.01*maximum(kappa(:)))
FileIO.save("kappa.jld2", "kappa", kappa)

#set up denominator
#denom = A' * (reshape(Base._sum(A', 1), size(wi)) .* wi)
thesum = A*(ones(420, 420, 96) .* mask)
denom = A' * (thesum .* wi)
FileIO.save("denom.jld2", "denom", denom)

function div0(a, b)
    if b != 0
        return a/b
    end
    return zero(typeof(b))
end
