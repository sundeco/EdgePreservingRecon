#combining matlab poisson.m, poisson1.m, poisson2.m
using Distributions
using SpecialFunctions

function poisson(xm ; factor = 0.85)
    dim = size(xm)
    xm = xm[:]
    try_poissrnd = false
    data = xm
    small = xm .< 12
    data[small] = poisson1(xm[small])
    data[.!small] = poisson2(xm[.!small], factor)
    return reshape(data, dim)
end

function poisson1(xmean)
    data = zeros(length(xmean))
    i_do = ones(length(xmean))
    ee = exp.(xmean)

    while any(i_do[:] .!= 0)
        i_do = ee .>= 1
        data[i_do] = 1 .+ data[i_do]
        ee = ee .* rand(length(xmean)) .* i_do
    end

    data = data .- 1
    return data
end

function poisson2(xm, factor)
    data = zeros(size(xm))

    if any(xm[:] .< 0)
        error("negative poisson means")
    end

    tmp = xm[xm .> 0]
    data[xm .> 0] = poisson2_positive(tmp[:], factor)
end

#assumes xm is a column and the output is a column
function poisson2_positive(xm, factor)
    sx = sqrt.(2*xm)
    lx = log.(xm)
    #gx = xm .* lx - SpecialFunctions.lgamma.(1 .+ xm)
    gx = 0.0 * xm
    for i = 1:length(xm)
        gx[i] = xm[i] * lx[i] - SpecialFunctions.logabsgamma(1 + xm[i])[1]
    end

    data = zeros(size(xm))
    id = [1:length(xm);]

    while any(id .!= 0)
        Tss = sx[id]
        Tll = lx[id]
        Tgg = gx[id]
        Txx = xm[id]

        yy = zeros(size(id))
        em = zeros(size(id))
        ib = ones(size(id)) .== ones(size(id))

        while any(ib .== true)
            yy[ib] = tan.(pi * rand(length(ib)))
            em[ib] = Tss[ib] .* yy[ib] + Txx[ib]
            #ib = find(em < 0) #this does not seem to have a good julia version
            ib = findall(x -> x < 0, em)
        end

        em = floor.(em)
        #tt = factor * (1 .+ yy.*yy) .* exp.(em .* Tll - SpecialFunctions.logabsgamma.(em .+ 1)[1] - Tgg)
        tt = 0.0 * em
        for i = 1:length(em)
            tt[i] = factor * (1 + yy[i]^2) * exp(em[i] * Tll[i] - SpecialFunctions.logabsgamma(em[i]+1)[1] - Tgg[i])
        end
        if any(tt .> 1)
            error("factor is too large")
        end

        ig = rand(length(id)) .< tt
        data[id[ig[:]]] = em[ig]
        id = id[.!ig[:]]
    end
    return data
end

test = 0
if test == 1
    xm = 10 .^(range(3, stop=6, length=100))
    out = poisson2(xm, 0.85)

end

function testpoisson(n)
    xm = 10 .^(range(3, stop=6, length=n))
    return poisson2(xm, 0.85)
end
