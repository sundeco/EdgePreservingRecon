using MIRTjim: jim
include("Reg1.jl")
include("subset_start.jl")
"""
penalized weighted least squares estimation / image reconstruction
using relaxed linearized augmented lagrangian method with
(optionally relaxed) ordered subsets. (relaxed OS-LALM)

See ?. 201? IEEE T-MI by Hung Nien & J A Fessler
"Relaxed linearized algorithms for faster X-ray CT image reconstruction"
http://dx.doi.org/10.1109/TMI.201?


cost(x) = (y-Ax)' W (y-Ax) / 2 + R(x)

in
	x	[np 1]		initial estimate
	Ab	[nd np]		Gblock object (needs abs(Ab) method)
				or sparse matrix (implies nsubset=1)
	yi	[nb na]		measurements (noisy sinogram data)
	R	penalty		object (see Reg1.jl), can be []

option
	niter			# of iterations (default: 1)
	wi	[nb na]		weighting sinogram (default: [] for uniform)
	pixmax	[1] or [2]	max pixel value, or [min max] (default [0 inf])
	denom	[np 1]		precomputed denominator
	aai	[nb na]		precomputed row sums of |Ab|
	relax0	[1] or [2]	relax0 or (relax0, relax_rate) (default 1)
	rho			AL penalty parameter (default: [] for decreasing rho)
	alpha			over-relaxation parameter (default: 1.999)
	userfun	@		user-defined function handle (see default below)
					taking arguments (x, userarg{:})
	chat

out
	xs	[np niter]	iterates

Translated from pwls_ep_os_rlalm_2d.m
Copyright 2022-5-31 Jason Hu and Jeff Fessler, University of Michigan
"""
function pwls_ep_os_rlalm_3d(x, A, yi, R::Reg1 ;
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
        mat"cg = ct_geom('ge2', 'down', 1);"

        mat"ig = image_geom('nx', 420, 'dx', 500/512, 'nz', 96, 'dz', 0.625, 'down', 1);"

        mat"ig.mask = ig.circ > 0;"
        mat"A = Gcone(cg, ig, 'type', 'sf2');"
        mat"Ab = Gblock(A, $nblock);"
        mat"Ab = block_op(Ab, 'ensure');"
        mat"nblock = block_op(Ab, 'n');"
        @mput wi
        @mput yi
        @mget nblock
    else
        Ab = A
    end
    @show(nblock)

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
        @time li = Ab[iblock] * MIRT.embed(x, mask)
    end
    odimsize = size(li)
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
        @time zeta = scale * Ab[iblock]' * reshape(resid, odimsize)
        zeta = zeta[mask]
        @show(sum(zeta))
    end

    g = rho2(1) * zeta
    h = denom .* x - zeta

    #some struct about info
    mask_roi = mask[:,:,17:80]
    xtrue_msk = xtrue[mask_roi]
    SqrtPixNum = sqrt(sum(mask_roi .> 0))

    for iter = 1:niter
        #xold = x
        relax = relax0 / (1 + relax_rate * (iter-1))

        if length(xtrue) != 0
            #compute info
        end

        #loop over subsets
        for iset = 1:nblock
            k = nblock*(iter-1) + iset

            num = rho2(k) * (denom .* x - h) .+ (1-rho2(k)) * g
            den = rho2(k) * denom

            @time num = num + R.cgrad(x)
            @time den = den + R.denom(x)

            x = x - relax * num ./ den
            x = max.(x, pixmin)
            x = min.(x, pixmax)
            y = MIRT.embed(x, R.mask)[:,:,:,1]
            @show(sqrt(sum((y[:,:,17:80]-xtrue).^2/sum(R.mask))))
            #x = x[:,:,1]

            iblock = starts[iset]
            ia = iblock:nblock:na

            if usemat
                @mput iblock
                @mput x
                mat"li = Ab{iblock} * x"
                @mget li
            else
                @time li = Ab[iblock] * MIRT.embed(x, mask)
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
                @time zeta = scale * Ab[iblock]' * reshape(resid, odimsize)
                zeta = zeta[mask]
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
        y = MIRT.embed(x, R.mask)[:,:,:,1]
        println(sqrt(sum((y[:,:,17:80]-xtrue).^2/sum(R.mask))))
    end
    #return xs
    return x
end
