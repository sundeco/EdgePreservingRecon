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
