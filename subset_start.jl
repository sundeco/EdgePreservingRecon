function subset_start(nsubset)
    if length(nsubset) == 1
        starts = 1 + bit_reverse(nsubset)
    end
    nsubset = length(starts)

    if any(sort(start[:])) .!= [1:nsubset]
        error("missing subset")
    end
end

function bit_reverse(mm)
    nn = convert(Int, 2^ceil(log2(mm)))
    ii = zeros(nn)
    for i = 0:nn-1
        ii[i] = Base.parse(Int, reverse(Base.bin(Unsigned(i), 1, false)); base=2)
    end
    ii = ii[ii .< mm]
    return ii
end
