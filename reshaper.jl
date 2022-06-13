"""
reshape function that is more flexible and allows for multiples

in
    x       [*dim (Ld)]
    dim     short row or column
            if dim is '2d' then y is 2d with first n-1 dims collapsed

out
    y       [dim (Ld)]

Translated from reshaper.m in MIRT
Copyright 2022-5-31, Jason Hu and Jeff Fessler, University of Michigan
"""
function reshaper(x, dim)
    dim_i = size(x)

    if isa(dim, String) && dim == "2d"
        return reshape(x, (prod(dim_i[1:end-1]), dim_i[end]))
    end

    dim_e = dim_i[2:end]
    b = Tuple(dim)
    dim_o = (b...,dim_e...)
    y = reshape(x, dim_o)
end

test = 0
if test == 1
    x = reshape([1:210;], (30,7))
    dim = [2,3,5]
    y = reshaper(x, dim)
    @show(y)
end
