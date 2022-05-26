include("Reg1.jl")
include("subset_start.jl")
function pwls_ep_os_rlalm_2d(x, A, yi, R ;
    niter::Int = 1,
    wi = [], #weighting matrix for sinogram
    pixmax = [0 100000000],
    denom = [], #precomputed denominator
    aai = [], #precomputed row sums of |Ab|
    relax0 = 1,
    rho = [], #AL penalty parameter
    alpha::Float32 = 1.999, #over relaxation parameter
    chat::Bool = false
    )
    nblock = 1
    isave = []
    pixmax = 100000000
    scale_nblock = true
    update_even_if_denom_0 = true
    xtrue = 0 * x
    mask = ones(size(xtrue))

    starts = subset_start(nblock)
    relax0 = relax0[1]
    if length(relax0) == 1
        relax_rate = 0
    elseif length(relax0) == 2
        relax_rate = relax0[2]
    else
        error("relax")
    end

    if length(pixmax) == 2
        pixmin = pixmax[1]
        pixmax = pixmax[2]
    elseif length(pixmax) == 1
        pixmin = 0
        pixmax = pixmax
    else
        error("pixmax")
    end

    if length(rho) == 0
        rho2(k) = pi/(alpha*k) * sqrt(1-(pi/(2*(alpha*k)))^2) * (k > 1) + (k == 1)
    else
        rho2(k) = rho
    end

    #initialization line 151 in matlab
    iblock = nblock
    ia = iblock:nblock:na
    li = Ab[iblock] * x
    li = reshape(li, (nb, length(ia)))
    resid = wi[:,ia] .* (li - yi[:,ia])
    if scale_nblock
        scale = nblock
    else
        scale = na / length(ia)
    end
    zeta = scale * Ab[iblock]' * resid[:]

    g = rho(1) * zeta
    h = denom .* x - zeta

    #some struct about info
    xtrue_msk = xtrue[mask]
    SqrtPixNum = sqrt(sum(mask .> 0))

    for iter = 1:niter
        xold = x
        relax = relax0 / (1 + relax_rate * (iter-1))

        if length(xtrue) != 0
            #compute info
        end

        #loop over subsets
        for iset = 1:nblock
            k = nblock*(iter-1) + iset

            num = rho2(k) * (denom .* x - h) .+ (1-rho2(k)) * g
            den = rho2(k) * denom
            if length(R) != 0
                num = num + R.cgrad(x)
                den = den + R.denom(x)
            end
            x = x - relax * num ./ den
            x = max(x, pixmin)
            x = min(x, pixmax)

            iblock = starts(iset)
            ia = iblock:nblock:na

            li = Ab[iblock] * x
            li = reshape(li, (nb, length(ia)))
            resid = wi[:,ia] .* (li - yi[:,ia])

            if scale_nblock
                scale = nblock
            else
                scale = na / length(ia)
            end

            zeta = scale * Ab[iblock]' * resid[:]
            g = (rho2(k) * (alpha * zeta + (1-alpha)*g) + g) / (rho2(k)+1)
            h = alpha * (denom .* x - zeta) + (1-alpha) * h
        end
    end
    if chat
        df = 0.5 * sum(wi[:] .* (Ab * x - yi[:])) .^2
        rp = R.penal(x)
    end
end
