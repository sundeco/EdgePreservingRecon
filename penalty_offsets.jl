"""
function offsets = penalty_offsets(offsets, isize)

determine the 'offsets' argument needed by regularizers

in
	offsets		empty or string or array
			'1d' [1] (1d default)
			'2d:hvd' [1 nx nx+1 nx-1] (2d default)
			'3d:hvu' [1 nx nx*ny] (3d default)
			'3d:26'	all 13 neighbors in 3d
	isize [N]	image size
out
	offsets [L]	penalty offsets, e.g., [1 nx nx+1 nx-1]

Translated from penalty_offsets.m in MIRT
Copyright 2022-5-29, Jason Hu and Jeff Fessler, University of Michigan
"""
function penalty_offsets(offsets, isize)
    if length(offsets) == 0
        if length(isize) == 1
            if isize[1] <= 1
                offsets = []
            else
                offsets = [1]
            end
        elseif length(isize) == 2
            if isize[2] == 1
                if isize[1] <= 1
                    offsets = []
                else
                    offsets = [1]
                end
            else
                nx = isize[1]
                offsets = [1, nx, nx+1, nx-1]
            end
        elseif length(isize) == 3
            nx = isize[1]
            ny = isize[2]
            offsets = [1, nx, nx*ny]
        else
            error("only 2D and 3D done")
        end
    elseif isa(offsets, Array{String}) || isa(offsets, String)
        if offsets == ["ident", "identity", "I"]
            offsets = [0]
        elseif offsets == "1d"
            if maximum(isize) <= 1
                offsets = []
            else
                offsets = [1]
            end
        elseif offsets == ["2d,hv", "2d:hv"]
            nx = isize[1]
            ny = isize[2]
            if nx == 1 && ny == 1
                offsets = []
            elseif nx == 1 || ny == 1
                offsets = [1]
            else
                offsets = [1, nx]
            end
        elseif offsets == ["2d,hvd", "2d:hvd"]
            nx = isize[1]
            ny = isize[2]
            if nx == 1 && ny == 1
                offsets = []
            elseif nx == 1 || ny == 1
                offsets = [1]
            else
                offsets = [1, nx, nx+1, nx-1]
            end
        elseif offsets == "3d:hvu"
            nx = isize[1]
            ny = isize[2]
            offsets = [1, nx, nx*ny]
            if nx == 1 || ny == 1
                error("not done")
            end
        elseif offsets == "3d:26"
            if length(isize) != 3
                error("expects 3d image")
            end
            nx = isize[1]
            ny = isize[2]
            nz = isize[3]
            tmp = [-1:1;] .+ ([-1:1;]*nx)'
            offsets = [1; nx .+ [0,1,-1]; nx*ny .+ tmp[:]]
        else
            error("Bad offsets string")
        end
    else
        error("offsets needs to be a string or empty")
    end
    return offsets
end
