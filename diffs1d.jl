#=
Left finite differences for a 1d array
Similar to diffl.jl in MIRT.jl, but works for arbitrary offsets

2022-5-29 Jason Hu and Jeff Fessler, University of Michigan
=#
export diffs1d_map
using LinearMapsAA: LinearMapAA, LinearMapAM, LinearMapAO

"""
Apply left finite difference operator to input array x and store the result in
pre allocated array g
Takes difference between two entries with difference in indices = offset
Only works on 1D arrays

Option:
flip: is -1 by default, set it to 1 to get the absolute value

Example: x = [4,5,3,12,7], offset = 2
Then output g = [0,0,-1,7,4]
"""
function diffs1d!(g, x, offset::Int ; flip = -1)
    if length(g) != length(x)
        error("not the same length")
    end
    if offset >= length(g)
        error("offset too big")
    end
    g[1:offset] .= 0
    for i = offset+1:length(g)
        #macro @inbounds
        @inbounds g[i] = x[i] + flip * x[i-offset]
    end
    return g
end

"""
Adjoint of the difference operator, g is input, output stored in z
Last offset elements are set equal to g
"""
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

"""
Create an operator for taking finite differences
in
N: length of the vector
offset: Difference in indices with which to take a difference

out
LinearMapAA object for computing finite differences
"""
function diffs1d_map(N, offset ; T::Type = Float32)
    forw!(g,x) = diffs1d!(g, x, offset)
    back!(z,g) = diffs1d_adj!(z, g, offset)

    return LinearMapAA(forw!, back!, (1,1) .* N)
end

"""
Create the absolute value of the operator for taking finite differences
Used in gradient and denominator of regularizer functions
in
N: length of the vector
offset: Difference in indices with which to take a sum

out
LinearMapAA object for computing finite sums
"""
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
