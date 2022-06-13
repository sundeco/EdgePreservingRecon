"""
Compute array of subset starting indices "starts" for OS algorithms.
If input is an empty matrix, then 1 subset is used.
If input is a scalar power of 2 != 1,
then the "bit-reversal ordering" is used.
If input is an array, then it is simply checked for completeness
and returned.
(nsubset can be 1 to #views)

See guan:94:apa doi 10.1088/0031-9155/39/11/013

Translated from subset_start.m in MIRT
Copyright 2022-5-31, Jason Hu and Jeff Fessler, University of Michigan
"""
function subset_start(nsubset)
    if length(nsubset) == 1
        starts = 1 .+ bit_reverse(nsubset)
    end
    starts = convert(Array{Int}, starts)
    nsubset = length(starts)

    if any(sort(starts[:]) .!= [1:nsubset;])
        print(sort(starts[:]))
        error("missing subset")
    end
    return starts
end

function bit_reverse(mm)
    nn = convert(Int, 2^ceil(log2(mm)))
    bits = convert(Int, ceil(log2(mm+1)))
    ii = zeros(nn)
    for i = 0:nn-1
        ii[i+1] = Base.parse(Int, reverse(Base.bin(Unsigned(i), bits , false)); base=2)
    end
    ii = ii[ii .< mm]
    return ii
end
