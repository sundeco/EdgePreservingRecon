using MIRT
include("penalty_offsets.jl")
include("potential_fun.jl")
include("Rweights.jl")
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
    if length(type_wt) == 0
        if length(kappa) <= 128^2
            type_wt = "pre"
        else
            type_wt = "fly"
        end
    end

    #line 122 in reg1.m not done yet
    offsets = penalty_offsets(offsets, size(kappa))
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
        push!(C1s, MIRT.diffl_map(dim, 1))
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
        wt = pot.wpod(d)
        wt = wt .* R.wt.col(mm)
        tmp = R.C1s[mm]' * (wt .* d)
        cgrad = cgrad + tmp
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
    return Reg1_mat_denom_sqs1_cell(R.C1s, R.mask, R.pot, R.wt, x)
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
    return penal
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
        cgrad[:,ll] = R.cgrad1_fun(R, x[:,ll])
    end
    return cgrad
end

Reg1_fun0 = Dict([
    (:C1, R -> Reg1_com_C1(R)),
    (:C, R -> Reg1_com_C(R)),
    (:penal, R -> (x -> Reg1_com_penal(R,x))),
    (:cgrad, R -> (x -> Reg1_com_cgrad(R,x))),
    (:denom_sqs1 , R -> (x -> Reg1_mat_denom_sqs1(R,x))),
    (:denom, R -> (x -> Reg1_mat_denom(R,x))),
    (:cgrad1_fun, R -> (x -> Reg1_cgrad1_fun(R,x))),
    (:dercurv, R -> (x -> Reg1_dercurv(R,x)))
])

Base.getproperty(R::Reg1, s::Symbol) =
    haskey(Reg1_fun0, s) ? Reg1_fun0[s](R) :
    getfield(R, s)

Base.propertynames(R::Reg1) =
    (fieldnames(typeof(R))..., keys(Reg1_fun0)...)
