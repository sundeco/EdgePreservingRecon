#Left finite differences only for 1D arrays of any offset
export diffs1d_map
using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO

function diffs1d!(g, x, offset ; flip = -1)
    if length(g) != length(x)
        error("not the same length")
    end
    if offset >= length(g)
        error("offset too big")
    end
    g[1:offset] .= 0
    for i = offset+1:length(g)
        g[i] = x[i] + flip * x[i-offset]
    end
    return g
end

function diffs1d_adj!(z, g, offset ; flip = -1)
    if length(z) != length(g)
        error("not same length")
    end
    if offset >= length(g)
        error("offset too big")
    end
    z[length(g)-offset+1:length(g)] = g[length(g)-offset+1:length(g)]
    for i = 1:length(g)-offset
        z[i] = g[i] + flip * g[i+offset]
    end
    return z
end

function diffs1d_map(N, offset ; T::Type = Float32)
    forw!(g,x) = diffs1d!(g, x, offset)
    back!(z,g) = diffs1d_adj!(z, g, offset)

    return LinearMapAA(forw!, back!, (1,1) .* N)
end

function diffs1d_map_abs(N, offset ; T::Type = Float32)
    forw!(g,x) = diffs1d!(g, x, offset ; flip = 1)
    back!(z,g) = diffs1d_adj!(z, g, offset ; flip = 1)

    return LinearMapAA(forw!, back!, (1,1) .* N)
end

test = 0
if test == 1
    x = [4,5,3,12,7]
    g = similar(x)
    z = similar(x)
    out = diffs1d!(g, x, 2)
    out2 = diffs1d_adj!(z, x, 2)

    A = diffs1d_map(5, 2)
    out = A*x
end
