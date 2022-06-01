using MIRT
include("penalty_offsets.jl")
include("potential_fun.jl")
include("Rweights.jl")
include("diffs1d.jl")
export Reg1

mutable struct Reg1
    type_penal::String
    edge_type::String
    type_diff
    pot_arg
    beta
    type_denom::String
    pre_denom_sqs1_x0::Bool
    distance_power::Float32
    user_wt
    order::Float32
    mask
    offsets
    offsets_is_zxy::Bool
    control::Float32
    type_wt::String
    nthread::Int
    M::Int
    dim
    isquad::Bool
    np::Int
    wt::Rweight
    C1s #should be some difference matrix thing
    pot #array of potential functions
    #cdp_arg #this appears to only be for mex, so ignore it
end

function Reg1(kappa ;
    type_penal::String = "",
    edge_type::String = "tight",
    type_diff::String = "",
    beta = 1,
    pot_arg = ["quad"], #huber, delta, cell array
    pot_arg2 = [1],
    type_denom::String = "none",
    pre_denom_sqs_x0::Bool = false,
    distance_power = 1,
    user_wt = [],
    order::Int = 1,
    mask = [],
    offsets = [],
    offsets_is_zxy::Bool = false,
    control = [],
    type_wt::String = "",
    nthread = 1)
#=
Build roughness penalty regularization "object" based on Cdiff1() objects,
for regularized solutions to inverse problems.
This version supercedes Robject() by providing its capabilities while
also providing options that use less memory.  By default it tries to use
a mex version (penalty_mex.mex*) but if that fails it reverts to a pure
matlab form that should be completely portable.

General form of (possibly nonquadratic) penalty function:
	R(x) = \sum_{m=1}^M sum_n w[n;m] potential_m( [C_m x]_n )
where M is the number of neighbors (offsets), and n=1,...,N=numel(x)
and C_m is a (N x N) differencing matrix in the direction off_m.

Penalty gradient is \sum_m C_m' D_m(x) C_m x,
where D_m(x) = diag{w[n;m] \wpot_[n;m]([C_m x]_n)} and \wpot(t) = \pot(t) / t

in
	kappa	[(N)]		kappa array, or logical support mask

options
	'type_penal'		'mat' only currently
	'type_diff'		'def|ind|mex|sparse|...' (see Cdiff1)
	'order', 1|2		1st-order only
	'offsets', [M] | char
				offsets to neighboring pixels
					(see Cdiff1 for the defaults)
				use '3d:26' to penalize all 13 pairs of nbrs
				use "0" for C = I (identity matrix)
	'beta', [1] | [M]	global regularization parameter(s)|				default: 2^0
	'pot_arg', {} 		arguments to potential_fun()
					e.g., {'huber', delta}, or cell{M} array
				default: {'quad'} for quadratic regularization.
 ?	'pre_denom_sqs1_x0'	precompute denominator for SQS at x=0? (def: 0)
	'type_denom', ''	type of "denominator"
					(for quadratic surrogates like SPS)
					todo: improve documentation!
		'matlab'	denominator for SPS
					todo: precompute?  or just on the fly?
		'aspire'	denominator for SPS that matches aspire
		'none'		no denominator precomputation (default)
	'distance_power', 0|1|2	See Rweights.jl
	'user_wt', [(N) M]	""
	'mask'			Override default: mask = (kappa ~= 0)
				(Only use if you sure know what you're doing!)

out
	Reg1  struct object with methods:
	R.penal(x)	evaluates R(x)
	R.cgrad(x)	evaluates \cgrad R(x) (column gradient)
	R.denom_sqs1(x)	evaluates denominator for separable quadratic surrogate
				\sum_m |C_m|' |D_m(x)| (|C_m| 1)
	R.denom(x)	evaluates denominator for separable surrogate
	[pderiv pcurv] = feval(R.dercurv, R, R.C1*x) derivatives and curvatures
				for non-separable parabola surrogates
	R.diag		diagonal of Hessian of R (at x=0), for preconditioners.
	R.C1		finite differencing matrix; for 1st-order differences
			nonzero entries of each row are 1 and -1;
			for 2nd-order differences each row contains {-1 2 -1}.
			Almost always should be used in conjunction with R.wt.
	R.C		diag(sqrt(w)) * C1, only for quadratic case!
	R.wt		see Rweights
=#
    if length(type_wt) == 0
        if length(kappa) <= 128^2
            type_wt = "pre"
        else
            type_wt = "fly"
        end
    end

    #line 122 in reg1.m not done yet
    offsets = penalty_offsets(offsets, size(kappa))
    offsets = convert(Array{Int}, offsets)

    M = length(offsets)

    if length(kappa) < order + 1
        error("image size too small")
    end

    if length(control) == 0
        if order == 0
            control = 0
        elseif order == 1
            control = 1
        elseif order == 2
            if length(user_wt) == 0
                control = 2
            else
                control = 1
            end
        else
            error("Order not done")
        end
    end

    dim = size(kappa)
    if ndims(kappa) == 2 && size(kappa,2) == 1 #1D case
        dim = size(kappa,1)
    end

    #line 165 of matlab file: R.kappa = kappa (why?)
    if !isreal(kappa)
        error("kappa must be real")
    end
    if any(kappa .< 0)
        error("kappa has negative values")
    end

    #Potential function setup, cell stuff with matlab
    isquad = true
    for mm = 1:M
        isquad = isquad && (pot_arg[min(mm, length(pot_arg))][1] == "quad")
    end

    if length(mask) == 0
        mask = (kappa .!= 0)
    else
        #check if mask is logical
    end
    np = sum(mask)

    wt = Rweights(kappa, offsets; type_wt, edge_type, beta, order, distance_power, user_wt)

    #differencing objects: look into diffl
    C1s = []
    for mm = 1:M
        #push!(C1s, MIRT.diffl_map(dim, 1))
        push!(C1s, diffs1d_map(prod(dim), offsets[mm]))
    end

    #desired potential function handles
    pot = []
    for mm = 1:M
        tmp = pot_arg[min(mm, length(pot_arg))]
        tmp2 = pot_arg2[min(mm, length(pot_arg2))]
        push!(pot, potential_fun(tmp, tmp2)) #this function is 974 lines lmao
    end
    type_penal = "mat"

    #line 234 of matlab code
    if type_penal == "mat"
        #don't use another function for this to avoid passing in 1000 arguments
        #R = Reg1_setup_mat()
        #original function starts line 274 of matlab
        if pre_denom_sqs_x0
            error("shouldn't run this)")
            #do something with denom_sqs_x0
        end
        R = Reg1(type_penal,
        edge_type, type_diff,
        pot_arg,
        beta,
        type_denom,
        pre_denom_sqs_x0,
        distance_power,
        user_wt,
        order,
        mask,
        offsets, offsets_is_zxy,
        control,
        type_wt,
        nthread,
        M,
        dim,
        isquad,
        np,
        wt,
        C1s, #should be some difference matrix thing
        pot)
    elseif type_penal == "mex"
        error("mex for julia omegalul")
    elseif type_penal == "zxy"
        error("not done yet")
    else
        error("bad penalty type")
    end
    return R
end

function Reg1_cgrad1_fun(R::Reg1, x)
    return Reg1_mat_cgrad1(R, x)
end

function Reg1_dercurv(R::Reg1, x)
    return 0
    #i don't know what's going on here but matlab doesn't have this function
end

function Reg1_mat_cgrad1(R::Reg1, x)
    cgrad = 0
    for mm = 1:R.M
        d = R.C1s[mm] * x
        pot = R.pot[mm]
        wt = pot.wpot(d)
        wt = wt .* R.wt.col(mm)
        tmp = R.C1s[mm]' * (wt .* d)
        cgrad = cgrad .+ tmp
    end
    cgrad = cgrad .* R.mask[:]
    return cgrad
end

function Reg1_mat_denom(R::Reg1, x)
    return Reg1_mat_denom_sqs1(R, x)
end

function Reg1_mat_denom_sqs1(R::Reg1, x)
    if length(x) == 0
        x = zeros(size(R.mask))
    end
    return Reg1_mat_denom_sqs1_cell(R.C1s, R.mask, R.pot, R.wt, x, R.offsets)
end

#located in a separate matlab file
function Reg1_mat_denom_sqs1_cell(C1s, mask, pots, ws, x, offsets)
    x = embed(x, mask)
    x = x[:]
    denom = 0
    for mm = 1:length(C1s)
        Cm = C1s[mm]
        d = Cm*x
        #Cm = abs(Cm) #change this
        Cm = diffs1d_map_abs(length(x), offsets[mm])
        ck = Cm * ones(size(x))
        pot = pots[mm]
        wt = pot.wpot(d)
        wt = wt .* reshape(ws.col(mm), size(wt))
        tmp = Cm' * (wt .* ck)
        denom = denom .+ tmp
    end
    denom = denom[mask[:]]
    return denom
end

function Reg1_com_penal(R::Reg1, x)
    if size(x,1) == R.np
        #matlab: x = embed(x, sr.mask), mask always has as many ones as the length of x
        x = MIRT.embed(x, R.mask)
    end
    x = reshape(x, (prod(R.dim), 1))
    LL = size(x, 2) #won't this always return 1?
    penal = zeros(LL, 1)
    for ll = 1:LL
        penal[ll] = Reg1_com_penal1(R, x[:,ll])
    end
    return penal[1]
end

function Reg1_com_penal1(R::Reg1, x)
    penal = 0
    for mm = 1:R.M
        d = R.C1s[mm] * x
        pot = R.pot[mm]
        d = pot.potk(d)
        wt = R.wt.col(mm)
        penal = penal + sum(wt .* d)
    end
    return penal
end

function Reg1_com_cgrad(R::Reg1, x)
    siz = size(x)
    x = MIRT.embed(x, R.mask)
    x = reshape(x, (prod(R.dim),1))
    LL = size(x, 2)

    cgrad = zeros(size(x))
    for ll = 1:LL
        cgrad[:,ll] = Reg1_cgrad1_fun(R, x[:,ll])
    end
    #ei.column does not work here, so just assume it's true
    cgrad = cgrad[R.mask[:], :]
    return cgrad
end

function Reg1_com_egrad(R::Reg1, x, delta)
    egrad = zeros(size(x))
    ej = zeros(size(x))
    for jj = 1:length(x)
        ej[jj] = 1
        egrad[jj] = (R.penal(x+delta*ej) - R.penal(x))/delta
        ej[jj] = 0
    end
    return egrad
end

Reg1_fun0 = Dict([
    (:C1, R -> Reg1_com_C1(R)),
    (:C, R -> Reg1_com_C(R)),
    (:penal, R -> (x -> Reg1_com_penal(R,x))),
    (:cgrad, R -> (x -> Reg1_com_cgrad(R,x))),
    (:egrad, R -> ((x,delta) -> Reg1_com_egrad(R,x, delta))),
    (:denom_sqs1 , R -> (x -> Reg1_mat_denom_sqs1(R,x))),
    (:denom, R -> (x -> Reg1_mat_denom(R,x))),
    #(:cgrad1_fun, R -> (x -> Reg1_cgrad1_fun(R,x))),
    (:dercurv, R -> (x -> Reg1_dercurv(R,x)))
])

Base.getproperty(R::Reg1, s::Symbol) =
    haskey(Reg1_fun0, s) ? Reg1_fun0[s](R) :
    getfield(R, s)

Base.propertynames(R::Reg1) =
    (fieldnames(typeof(R))..., keys(Reg1_fun0)...)
