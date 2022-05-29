using MIRTjim: jim
include("Reg1.jl")
include("subset_start.jl")
function pwls_ep_os_rlalm_2d(x, A, yi, R::Reg1 ;
    niter::Int = 1,
    wi = [], #weighting matrix for sinogram
    pixmax = Inf,
    denom = [], #precomputed denominator
    aai = [], #precomputed row sums of |Ab|
    relax0 = 1,
    rho = [], #AL penalty parameter
    alpha = 1.999, #over relaxation parameter
    chat::Bool = false,
    usemat::Bool = false,
    isave::String = "",
    nblock::Int = 1,
    xtrue = 0 * x,
    mask = ones(size(xtrue))
    )
    scale_nblock = true
    update_even_if_denom_0 = true

    if isave == "last"
        isave = niter
    else
        isave = 0:niter
    end
    #all of this corresponds to line 71 in matlab file
    if usemat
        if !@isdefined(irtdir)
        	ENV["MATLAB_ROOT"] = "/Applications/matlab"

        	irtdir = "/Users/jasonhu/Documents/MATLAB/phd/irt"
        	tmp = "addpath('$irtdir')"
        	eval_string(tmp)
        	mat"setup"
        end
        mat"sg = sino_geom('ge1', 'units', 'mm', 'strip_width', 'd', 'down', 1);"

        mat"ig = image_geom('nx', 420, 'dx', 500/512, 'down', 1);"

        mat"ig.mask = ig.circ > 0;"
        mat"A = Gtomo2_dscmex(sg, ig);"
        mat"Ab = Gblock(A, $nblock);"
        mat"nblock = block_op(Ab, 'n');"
        @mput wi
        @mput yi
        @mget nblock
    end

    starts = subset_start(nblock)

    if length(wi) == 0
        wi = ones(size(yi))
    end
    if length(aai) == 0 && usemat
        mat"aai = reshape(sum(abs(Ab)'), size(yi));"
        @mget aai
    end

    if ndims(yi) != 2 || (size(yi, 2) == 1 && nblock > 1)
        error("bad yi size")
    end
    if ndims(wi) != 2 || (size(wi, 2) == 1 && nblock > 1)
        error("bad wi size")
    end

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

    #line 114 of matlab
    if length(denom) == 0 && usemat
        mat"denom = abs(Ab)' * col(aai .* wi);"
        @mget denom
    end

    if alpha < 1 || alpha > 2
        error("alpha should be between 1 and 2")
    end
#=
    if length(rho) == 0
        rho2(k) = pi/(alpha*k) * sqrt(1-(pi/(2*(alpha*k)))^2) * (k > 1) + (k == 1)
    else
        rho2(k) = rho
    end
=#
    rho2(k) = pi/(alpha*k) * sqrt(1-(pi/(2*(alpha*k)))^2) * (k > 1) + (k == 1)

    (nb,na) = size(yi)
    x = x[:]
    np = length(x)
    xs = zeros(np, length(isave))
    if any(isave .== 0)
        xs[:, isave .== 0] = x
    end

    #initialization line 151 in matlab
    iblock = nblock
    ia = iblock:nblock:na
    if usemat
        @mput iblock
        @mput x
        mat"li = Ab{iblock} * x"
        @mget li
    else
        li = Ab[iblock] * x
    end
    li = reshape(li, (nb, length(ia)))
    resid = wi[:,ia] .* (li - yi[:,ia])
    if scale_nblock
        scale = nblock
    else
        scale = na / length(ia)
    end
    if usemat
        @mput resid
        @mput iblock
        @mput scale
        mat"zeta = scale * Ab{iblock}' * resid(:)"
        @mget zeta
    else
        zeta = scale * Ab[iblock]' * resid[:]
    end

    g = rho2(1) * zeta
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

            num = num + R.cgrad(x)
            den = den + R.denom(x)

            x = x - relax * num ./ den
            x = max.(x, pixmin)
            x = min.(x, pixmax)
            x = x[:,:,1]

            iblock = starts[iset]
            ia = iblock:nblock:na

            if usemat
                @mput iblock
                @mput x
                mat"li = Ab{iblock} * x"
                @mget li
            else
                li = Ab[iblock] * x
            end
            li = reshape(li, (nb, length(ia)))
            resid = wi[:,ia] .* (li - yi[:,ia])

            if scale_nblock
                scale = nblock
            else
                scale = na / length(ia)
            end

            if usemat
                @mput resid
                @mput iblock
                @mput scale
                mat"zeta = scale * Ab{iblock}' * resid(:)"
                @mget zeta
            else
                zeta = scale * Ab[iblock]' * resid[:]
            end
            g = (rho2(k) * (alpha * zeta + (1-alpha)*g) + g) / (rho2(k)+1)
            h = alpha * (denom .* x - zeta) + (1-alpha) * h
        end

        # if any(isave .== iter)
        #     xs[:, isave .== iter] = x
        # end
        if chat
            if usemat
                @mput wi
                @mput x
                @mput yi
                df = 0.5 * sum(wi(:) .* (Ab * x - yi(:))) .^2
                @mget df
            else
                df = 0.5 * sum(wi[:] .* (Ab * x - yi[:])) .^2
            end
            rp = R.penal(x)
        end
        println(iter)
        y = MIRT.embed(x, R.mask)[:,:,1]
        println(assess_ssim(y, xtrue))
    end
    #return xs
    return x
end
