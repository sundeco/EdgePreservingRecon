using MIRT
include("penalty_displace.jl")
function penalty_distance(offsets, sizes)
#=
convert scalar offsets to Euclidean distances to neighbors,
with the exception that '0' is mapped to 1 for identity case

 in
	offsets		[MM]
	sizes		[1 ndim]
out
	distance	[MM 1]

Translated from penalty_distance.m
Copyright 2022-5-23, Jason Hu and Jeff Fessler, University of Michigan
=#

    displace = penalty_displace(offsets[:], sizes)
    distance = sqrt.(sum(displace.^2, dims = 2))
    distance[offsets .== 0] .= 1
    return distance
end

test = 0
if test == 1
    nx = 100
    ny = 80
    ix, iy = ndgrid(-2:2, -2:2)
    tmp = ix + iy * nx
    offsets = tmp[:]
    dis = penalty_distance(offsets, [nx, ny])
    tru = sqrt.(ix[:].^2 + iy[:].^2)
    #tru[(ix .== 0) .* (iy .== 0)] .= 1 i'm too lazy to debug this
    tru[13] = 1
    print(dis-tru)
end
