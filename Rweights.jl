include("penalty_displace.jl")
include("penalty_distance.jl")
using MIRT
using MIRTjim
using MAT #only for testing purposes
using ThreadsX

export Rweight

mutable struct Rweight
    offsets
    type_wt::String
    beta
    edge_type::String
    distance_power::Float32
    order::Float32
    user_wt
    isize
    beta_dis
    displace
    kappa
    funtype::String #used to store what function col will call
end

"""
Build "weights" for roughness penalty for regularized methods.
Intended to be used internally by roughness penalty object Reg1.

General form of roughness penalty:
	R(x) = sum_{m=1}^M sum_n w[n;m] potential( [C_m x]_n )
where M is the number of neighbors (offsets), and n=1,...,N=numel(x)
and C_m is a (square) differencing matrix in the mth direction off_m.

General form of weights:
	w[n;m] = beta_m / || off_m ||^distance_power kappa2[n;m] user_wt[n;m]
where form of kappa2[n;m] depends on 'edge_type' as follows:
	kappa^2[n] for 'simple' case,
	kappa[n] kappa[n-off_m] for 'tight', order = 1
	kappa[n] sqrt(kappa[n-off_m] * kappa[n+off_m]) for 'tight', order = 2
	kappa[n] kappa[n-off_m] for 'leak' case, order = 1,
			unless either is zero, in which case square the other.
			for order=2, square the maximum of all three if needed.
	kappa[n] kappa_extend[n-off_m] for 'aspire2' case, order = 1

Although the 'tight' case seems preferable to avoid any influence
of pixels outside of the support, the 'simple' case has the advantage
that it can be precomputed easily with only N storage, instead of M*N.

in
	kappa	[(N)]		kappa array, or logical support mask
	offsets	[M]		offset(s) to neighboring pixels

options
	'type_wt'		what type of object to return:
		'array'		[(N) M] array of w[n;m] values
		'fly'		strum object for computing w[:,m] on-the-fly
		'pre'		strum object for precomputed w[:,m] (default)
	'edge_type'		how to handle mask edge conditions
		'simple'	kappa^2 (almost "tight" but saves memory)
		'tight'		only penalize within-mask differences (default)
		'leak'		penalize between mask and neighbors
					(mostly for consistency with ASPIRE)
	'beta', [1] or [M]	regularization parameter(s) (default: 2^0)
	'order' {1|2}		differencing order (only for some cases)
	'distance_power', {0:1:2}  1 classical (default), 2 possibly improved
					use 0 if user_wt has its effect already
	'user_wt', [(N) M]	User-provided array of penalty weight values
				of dimension [size(kappa) length(offsets)].

out
	wt [(N) M]		w[n;m] values needed in regularizer
	or a Rweight object with the following public methods:
		wt.col(m)	returns w[:,m]

Translated from Rweights.m in MIRT
Copyright 2022-5-31, Jason Hu and Jeff Fessler, University of Michigan
"""
function Rweights(kappa, offsets ;
    type_wt::String = "pre", #should be either array, fly, or pre
    edge_type::String = "tight",
    beta = 1, #[1] or [M]
    order::Int = 1, #differencing order 1 or 2
    distance_power::Int = 1, #0 1 or 2
    user_wt = [],
    )
    isize = size(kappa)

    #beta values adjusted for distance
    if length(beta) == 1
        beta_dis = beta ./ penalty_distance(offsets[:], isize) .^ distance_power
    else
        beta_dis = beta[:] ./ penalty_distance(offsets[:], isize) .^ distance_power
    end

    displace = penalty_displace(offsets, isize)

    #skip the aspire2 testing option
    if type_wt == "array"
        wt = Rweights_array(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
    elseif type_wt == "fly"
        wt = Rweights_strum_fly(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
    elseif type_wt == "pre"
        wt = Rweights_strum_pre(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
    else
        error("unknown type")
    end
    return wt
end

function Rweights_array(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
    MM = length(offsets)
    #line 221 of matlab file
    wt = zeros(prod(isize), MM)
    for mm = 1:MM
        wt[:,mm] = beta_dis(mm) * Rweights_kappa2_mat(kappa, offsets[mm], displace[mm,:], edge_type, order)
    end
    #incorporate user provided weights
    if length(user_wt) != 0
        tmp = reshape(user_wt, size(wt))
        wt = wt .* tmp
    end
    return wt
end

function Rweights_strum_pre(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
    return Rweight(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, "column of w")
end

function Rweights_strum_fly(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
    if edge_type == "simple"
        #kappa2 = kappa[:].^2 don't need this here because it's automatically included
        funtype = "strum fly col simple"
    else
        funtype = "strum fly col array"
    end
    if length(user_wt) != 0
        funtype = "strum mul user wt"
    end
    return Rweight(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, funtype)
end

function Rweights_strum_fly_col_array(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, mm)
    return beta_dis[mm] * Rweights_kappa2_mat(kappa, offsets[mm], displace[mm, :], edge_type, order)
end

function Rweights_strum_mul_user_wt(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, mm)
    error("not done yet")
end

function Rweights_strum_fly_col_simple(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, mm)
    return beta_dis[mm] * kappa[:].^2
end

function Rweights_strum_show(wt::Rweight)
    #error("not done yet")
    isize = wt.isize
    offsets = wt.offsets
    Nd = isize
    Ns = prod(Nd)
    MM = length(offsets)
    tmp = zeros(Ns, MM)
    for mm = 1:MM
        tmp[:,mm] = wt.col(mm)
    end

    tmp = reshape(tmp[:,1], Nd)
    jim(tmp)
end

function Rweights_column(wt::Rweight, l::Int)
    offsets = wt.offsets
    type_wt = wt.type_wt
    beta = wt.beta
    edge_type = wt.edge_type
    distance_power = wt.distance_power
    order = wt.order
    user_wt = wt.user_wt
    isize = wt.isize
    beta_dis = wt.beta_dis
    displace = wt.displace
    kappa = wt.kappa
    if wt.funtype == "column of w"
        w = Rweights_array(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa)
        return w[:,l]
    elseif wt.funtype == "strum fly col simple"
        return Rweights_strum_fly_col_simple(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, l)
    elseif wt.funtype == "strum fly col array"
        return Rweights_strum_fly_col_array(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, l)
    elseif wt.funtype == "strum mul user wt"
        return Rweights_strum_mul_user_wt(offsets, type_wt, beta, edge_type, distance_power, order, user_wt, isize, beta_dis, displace, kappa, l)
    else
        error("invalid function type")
    end
end

function Rweights_kappa2_mat(kappa, offset, displace, edge_type, order)
    Ns = length(kappa)
    Nd = size(kappa)
    displace = reshape(displace, (1, length(displace)))
    #pn = jf_protected_names

    if edge_type == "simple"
        kappa2 = kappa[:].^2
    elseif edge_type == "leak"
        kappa2 = zeros(size(kappa[:]))
        if order == 1
            ii = (1 + max(offset,0)):(Ns + min(offset,0))
            kappa0 = kappa[ii] + kappa[ii .- offset] .* (1 .- kappa[ii])
            kappa1 = kappa[ii .- offset] + kappa[ii] .* (1 .- kappa[ii .- offset])
            kappa2[ii] = kappa0 .* kappa1

            #something about pn and edge effects?
            ic = jf_ind2sub(Nd, ii)
            inn = jf_ind2sub(Nd, ii .- offset)
            good = (ic-inn) .== repeat(displace, inner = (1,1), outer = (length(ii), 1))
            #in matlab, all(good, 2) checks each row for being all nonzero
            kappa2[ii] = kappa2[ii] .* all(good)
        elseif order == 2
            ii = [(1 .+ abs.(offset)) : (Ns .- abs.(offset))]
            kaps = [kappa[ii] kappa[ii .- offset] kappa[ii .+ offset]]
            tmp = kaps[:,1] .* sqrt.(kaps[:,2] .* kaps[:,3])
            bad = (tmp .== 0)
            tmp[bad] = maximum(kaps[bad,:], dims=2) .^2
            kappa2[ii] = tmp

            #handle edge effects
            ic = jf_ind2sub(Nd, ii)
            inn = jf_ind2sub(Nd, ii .- offset)
            ip = jf_ind2sub(Nd, ii .+ offset)
            tmp = repeat(displace, inner = (1,1), outer = (length(ii), 1))
            good = ((ic-inn) .== tmp .&& (ip-ic) .== tmp)
            @inbounds kappa2[ii] = kappa2[ii] .* all(good, dims=2)
        else
            error("bad order")
        end
    elseif edge_type == "tight"
        kappa2 = zeros(size(kappa[:]))
        if order == 1
            ii = (1 + max(offset,0)):(Ns + min(offset,0)) #length is around 420*420*96 usually
            @inbounds kappa2[ii] = kappa[ii] .* kappa[ii .- offset]
            #ic = jf_ind2sub2(Nd, ii) #array of size 1.6e7 rows and 3 columns, there must be a faster way!
            #inn = jf_ind2sub2(Nd, ii .- offset)
            ic = jason_ind2sub(Nd[1], Nd[2], Nd[3], ii[1], ii[end])
            inn = jason_ind2sub(Nd[1], Nd[2], Nd[3], ii[1]-offset, ii[end]-offset)
            #good = (ic - inn) .== repeat(displace, inner = (1,1), outer = (length(ii), 1))

            #@time @inbounds good = [(@view(ic[i,:]) .== (@view(inn[i,:]) + displace)) for i = 1:size(ic,1)]
            #v = @view kappa2[ii]
            # for i = 1:size(ic,1)
            #     if all(j -> ic[i,j] == inn[i,j] + displace[j], 1:length(displace))
            #         v[i] = zero(eltype(kappa2))
            #     end
            # end
            # for i = 1:size(ic, 1)
            #     if !all(good[i])
            #         v[i] = zero(eltype(kappa2))
            #     end
            # end
            #good = dont_repeat(ic, inn, displace)
            #@inbounds kappa2[ii] .*= all(good, dims=2)
            get_kappa2!(ic, inn, displace, ii, kappa2)
        elseif order == 2
            ii = (1 + abs(offset)) : (Ns - abs(offset))
            kappa2[ii] = kappa[ii] .* sqrt.(kappa[ii .- offset] .* kappa[ii .+ offset])

            #handle edge effects
            ic = jf_ind2sub(Nd, ii)
            inn = jf_ind2sub(Nd, ii .- offset)
            ip = jf_ind2sub(Nd, ii .+ offset)
            tmp = repeat(displace, inner = (1,1), outer = (length(ii), 1))
            good = ((ic-inn) .== tmp .&& (ip-ic) .== tmp)
            kappa2[ii] = kappa2[ii] .* all(good, dims=2)
        else
            error("bad order")
        end
    elseif edge_type == "aspire2"
        if order != 1 || ndims(kappa) != 2
            error("bug YEP")
        end
        kappa = Rweights_kappa_expand(kappa)
        kappa2 = Rweights_kappa2_mat(kappa, offset, displace, "tight", order)
    else
        error("unknown edge type")
    end
    return kappa2
end

function Rweights_kappa_expand(kappa)
    tmp = kappa[[1, end], :]
    tmp = tmp[:]
    tmp2 = kappa[:,[1, end]]
    tmp2 = tmp2[:]
    if any(tmp) == 1 || any(tmp2) == 1
        error("kappa expand requires support one pixel from edge")
    end
    mask = kappa .!= 0
    (nx, ny) = size(mask)
    for iy = 2:ny-1
        for ix = 2:nx-1
            if mask[ix,iy] == 0
                continue
            end
            for dy = -1:1
                for dx = -1:1
                    if mask[ix+dx, iy+dy] == 0
                        kappa[ix+dx, iy+dy] = kappa[ix, iy]
                    end
                end
            end
        end
    end
    return kappa
end

function jf_ind2sub(Nd, ii)
    i2s = CartesianIndices(zeros(Nd))
    @inbounds out = i2s[ii]
    if length(Nd) == 2
        return [getindex.(out, 1) getindex.(out, 2)]
    elseif length(Nd) == 3
        return [getindex.(out, 1) getindex.(out, 2) getindex.(out, 3)]
    else
        error("not done")
    end
end

#get kappa2 directly
function get_kappa2!(ic, inn, displace, ii, kappa2)
    for i = 1:size(ic, 1)
        if ic[i, 1] - inn[i,1] != displace[1] || ic[i, 2] - inn[i,2] != displace[2]  || ic[i, 3] - inn[i,3] != displace[3]
            kappa2[ii[i]] = 0
        end
    end
    return kappa2
end

#replaces a line that requires repeating a matrix lots of times
function dont_repeat(ic, inn, displace)
    out = Vector{Int64}(undef, size(ic, 1))
    for i = 1:size(ic, 1)
        if ic[i, 1] - inn[i,1] == displace[1] && ic[i, 2] - inn[i,2] == displace[2]  && ic[i, 3] - inn[i,3] == displace[3]
            out[i] = 1
        else
            out[i] = 0
        end
    end
    return out
end

#version of jf_ind2sub2 that is faster for ranges of integers
function jason_ind2sub(m, n, p, lb, ub)
    first = jf_ind2sub2((m,n,p), [lb])
    out = Matrix{Int64}(undef, ub-lb+1, 3)
    out[1,1] = first[1]
    out[1,2] = first[2]
    out[1,3] = first[3]
    for i = 1:ub-lb
        if out[i,1] < m
            out[i+1, 1] = out[i,1] + 1
            out[i+1, 2] = out[i,2]
            out[i+1, 3] = out[i,3]
        elseif out[i,2] < n
            out[i+1,1] = 1
            out[i+1, 2] = out[i, 2] + 1
            out[i+1, 3] = out[i,3]
        else
            out[i+1,1] = 1
            out[i+1, 2] = 1
            out[i+1,3] = out[i,3] + 1
        end
    end
    return out
end

function jf_ind2sub2(Nd, ind)
    ind = ind[:]
    # subs = zeros(length(ind), length(Nd))
    if length(Nd) == 2
        first, second = Base._ind2sub(Nd, ind)
        # subs[:,1] = first
        # subs[:,2] = second
        return [first second]
    elseif length(Nd) == 3
        first, second, third = Base._ind2sub(Nd, ind)
        # subs[:,1] = first
        # subs[:,2] = second
        # subs[:,3] = third
        return [first second third]
    elseif length(Nd) == 4
        first, second, third, fourth = Base._ind2sub(Nd, ind)
        # subs[:,1] = first
        # subs[:,2] = second
        # subs[:,3] = third
        # subs[:,4] = fourth
        return [first second third fourth]
    else
        fail("not done")
    end
    return 0
end

function Rweights_self_test()

    return
end

Rweight_fun0 = Dict([
    (:show, wt -> Rweights_strum_show(wt)),
    (:col, wt -> (l::Int -> Rweights_column(wt, l)))
])

Base.getproperty(wt::Rweight, s::Symbol) =
    haskey(Rweight_fun0, s) ? Rweight_fun0[s](wt) :
    getfield(wt, s)

Base.propertynames(wt::Rweight) =
    (fieldnames(typeof(wt))..., keys(Rweight_fun0)...)

test = 0
if test == 1
    vars = matread("/Users/jasonhu/Documents/julia_files/edge_preserve/testRweights.mat")
    kappa = vars["kappa"]
    offsets = [1, 420, 421, 419]
    type_wt = "fly"
    edge_type = "tight"
    beta = 65536
    order = 1
    distance_power = 1
    user_wt = []

    mw = Rweights(kappa, offsets ;
        type_wt = type_wt,
        edge_type = edge_type,
        beta = beta,
        order = order,
        distance_power = distance_power,
        user_wt = user_wt)
    println(sum(mw.col(1)))
end
