using MIRT
"""
Convert scalar offsets to vector displacements, i.e.,
find d1,d2,... such that d1*1 + d2*n1 + d3*n1*n2 + ... = offset.
Example: if offset = n1, then [d1 d2 d3] = [0 1 0]
Assumes that d <= floor(n/2) for n > 2

in
	offsets	[LL]		see penalty_offsets.m
	sizes	[1 ndim]

out
	displace [LL ndim]	[d1 d2 ...] for each offset

Translated from penalty_displace.m in MIRT
Copyright 2022-5-23, Jason Hu and Jeff Fessler, University of Michigan
"""
function penalty_displace(offsets, sizes)
    ndim = length(sizes)
    offsets = offsets[:]
    displace = 6969 * ones(length(offsets), ndim)

    if length(sizes) == 2 && sizes[1] == 2 #special case
        if offsets[1] == 1 && offsets[2] == 2 && offsets[3] == 3 && offsets[4] == 1
            displace = [1 0; 0 1; 1 1; -1 1]
        elseif offsets[1] == 1 && length(offsets) == 1
            displace = [1,0]
        elseif offsets[1] == 2 && length(offsets) == 1
            displace = [0,1]
        elseif offsets[1] == 1 && offsets[2] == 2 && length(offsets) == 2
            displace = [1 0; 0 1]
        else
            error("not done")
        end
    elseif any(sizes[1:end-1] .== 2)
        error("not done: ambiguous")
    else
        residual = offsets
        for id = ndim:-1:1
            nd = prod(sizes[1:id-1])
            displace[:,id] = round.(residual ./ nd)
            residual -= displace[:,id] * nd
        end
        if any(residual .!= 0)
            error("bug kekw")
        end
    end

    #offset_check = displace * [1, cumprod(sizes[1:end-1])]
    #print(offsets == offset_check)
    return displace
end

test = 0
if test == 1
    nx = 2
    ny = 2
    offsets = [1, nx, nx+1, nx-1]
    dd = penalty_displace(offsets, [nx, ny])
    print(dd)

    nx = 2
    ny = 3
    dd = penalty_displace(offsets, [nx, ny])
    print(dd)

    nx = 3
    ny = 3
    offsets = [1, nx, nx+1, nx-1]
    dd = penalty_displace(offsets, [nx, ny])
    print(dd)

    nx = 100
    ny = 80
    ix, iy = ndgrid(-2:2, -2:2)
    tmp = ix + iy * nx
    offsets = tmp[:]
    dd = penalty_displace(offsets, [nx, ny])
    print(dd)
end
