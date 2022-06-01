export potential

mutable struct potential
    ptype::String
    delta
    param
end

function potential_fun(ptype::String, delta ; param = [])
#=
Define roughness penalty potential functions.
%|
The penalty will have the form
	R(x) = sum_k w_k * potential_k([Cx]_k, delta_k)
where w_k is provided elsewhere, not here!
%
in
	ptype		quad broken huber hyper2 hyper3 cauchy qgg2 genhub
			lange1 lange3 (fair) geman&mcclure gf1 li98cfs
			Recommended: 'hyper3'
			To list all possible choices, use:
			potential_fun('list') for "smooth" options
			potential_fun('list1') for "non-smooth" options
	delta		scalar, or image-sized array;
			"cutoff" parameter for edge-preserving regularization
	param		optional additional parameter(s) for some choices:
				'gf1' (generalized Fair) : [a b]
				'qgg2' : q
				'genhub' & 'stevenson94dpr' : [p q]
				'table1' : [dz, dpot([1:K] * dz)]

out
	pot		potential object, with data: ptype, delta and param
	methods:
		pot.potk(C*x)	potential function value
		pot.wpot(C*x)	potential 'weights' (aka half-quad. curvatures)
		pot.dpot(C*x)	potential derivative
 		pot.shrink(b, reg)	proximal (shrinkage) operation:
					argmin_z 1/2 |z - b|^2 + reg * pot(z)
		pot.plot()	plot all of the above functions

trick: for ptype 'gf1-fit', the param argument should be:
	{'name of potential to fit', points, param}
and this function returns the [a b] parameters needed for a subsequent
call with ptype 'gf1'

Translated from potential_fun.m in MIRT
Copyright 2022-5-29, Jason Hu and Jeff Fessler, University of Michigan
=#
    if ptype == "list1"
        return ["l0", "l1", "lpp", "tav", "broken", "fair-l1"]
    end
    scale = 1
    if ptype == "gf1-fit"
        return ir_potential_fun_gf1_fit(param)
    end
    if ptype == "hyper3"
        ptype == "hyper2"
        delta = delta/sqrt(3)
    end
    if ptype == "gf1" && param[1] == 0
        if param[2] != 1
            error("only b=1 makes sense for gf1 with a=0")
        end
        ptype = "lange3"
        param = []
    end
    return potential(ptype, delta, param)
end

function ir_potential_fun_gf1_fit(ptype, sv, param)
    sv = sv[:]
    pot = potential_fun(ptype, 1, param)
    pt = pot.wpot(sv)

    s1 = sv[1]
    s2 = sv[2]
    w1 = pt[1]
    w2 = pt[2]
    abs = [0,0]
    ab[1] = (w2 * s2 * (1-w1) -w1*s1*(1-w2)) / ((w1-w2)*s1*s2)
    ab[2] = (s2 * (1-w1) -s1*(1-w2)) / ((w1-w2)*s1*s2)
    return ab
end

# function ir_potential_fun_parse(ptype, delta, param)
#     #line 122 of matlab
#     dpot = []
#     shrink = []
#
#     if ptype == "lange3"
#
#     end
#     return ptype, delta, param, potk, wpot, dpot, shrink
# end

function potk(pot::potential, z)
    if pot.ptype == "lange3"
        return pot.delta.^2 .* (abs.(z ./ pot.delta) - log.(1 .+ abs.(z ./ pot.delta)))
    elseif pot.ptype == "quad"
        return (abs.(z).^2)/2
    elseif pot.ptype == "broken"
        return min.(z.^2, pot.delta.^2)/2
    elseif pot.ptype == "l0"
        return z .!= 0
    elseif pot.ptype == "abs" || pot.ptype == "l1"
        return abs.(z)
    else
        error("not done yet")
    end
end

function wpot(pot::potential, z)
    if pot.ptype == "lange3"
        return 1 ./ (1 .+ abs.(z ./ pot.delta))
    elseif pot.ptype == "quad"
        return ones(size(z))
    elseif pot.ptype == "broken"
        return abs.(z) .< pot.delta
    elseif pot.ptype == "l0"
        return nan(size(z))
    elseif pot.ptype == "l1" || pot.ptype == "abs"
        return 1 ./ abs.(z)
    else
        error("not done yet")
    end
end

function ir_potential_fun_shrink(pot::potential, b, reg)
    if pot.ptype == "lange3"
        out = zeros(size(b))
        for ii = 1:length(b)
            a = convert(Float32, abs(b[ii]))
            s = sign(b[ii])
            cost(z) = 0.5*(z[1]-a)^2 + reg * pot.potk(z[1])
            z0 = optimize(cost, [0.0])
            z0 = s*Optim.minimizer(z0)
            z1 = optimize(cost, [a])
            z1 = s*Optim.minimizer(z1)
            if cost(z0) < cost(z1)
                out[ii] = z0[1]
            else
                out[ii] = z1[1]
            end
        end
        return out
    elseif pot.ptype == "quad"
        return b ./ (1 .+ reg)
    elseif pot.ptype == "broken"
        return ir_broken_shrink(b, reg, pot.delta)
    elseif pot.ptype == "l0"
        return b .* (abs.(b) .> sqrt.(2*reg))
    elseif pot.ptype == "l1" || pot.ptype == "abs"
        return sign.(b) .* max.(abs.(b)-reg, 0)
    else
        error("not done yet")
    end
end

function dpot(pot::potential, z)
    if pot.ptype == "lange3"
        return z ./ (1 .+ abs.(z ./ pot.delta))
    elseif pot.ptype == "quad"
        return z
    elseif pot.ptype == "broken"
        return z .* (abs.(z) .< pot.delta)
    elseif pot.ptype == "l0"
        return nan(size(z))
    elseif pot.ptype == "abs" || pot.ptype == "l1"
        return sign.(z)
    else
        error("not done yet")
    end
end

function ir_potential_fun_plot(pot::potential)
    error("not done yet")
    z = LinRange(-1, 1, 101) * 2
    z = z * pot.delta
end

#this is probably unnecessary for julia
# function ir_tonumeric(x,y)
#     if isa(y, Float32)
# end

function ir_broken_shrink(z, reg, delta)
    out = z ./ (1 .+ reg)
    big = delta * (1 .+ reg) .< abs.(z)
    out[big] = z[big]
end

function ir_lpp_1_2_shrink(b, reg)
    sb = sign.(b)
    a = abs.(b)
    out = zeros(size(b))
    t0 = 2*sqrt.(a/3) .* cos.(1/3 * acos.(3*(reg/2)/2 ./ (-a) .* sqrt.(3 ./ a)))
    a_1_2 = 3/2 * reg ^ (2/3)
    big = a .> a_1_2
    out[big] = t0[big] .^ 2
    return out .* sb
end

function ir_lpp_4_3_shrink(b, reg)
    sb = sign.(b)
    a = abs.(b)
    x = sqrt.(a.^2 .+ 256 * reg^3 / 729)
    out = sb .* (a + 4*reg/(3*2^1/3) * ((x-a) .^ (1/3) - (x+a) .^ (1/3)))
    return out
end

function ir_huber_shrink(z, reg, delta)
    out = z ./ (1 .+ reg)
    big = delta .* (1 .+ reg) .< abs.(z)
    if length(reg) > 1
        reg = reg[big]
    end
    if length(delta) > 1
        delta = delta[big]
    end
    out[big] = z .* (1 .- reg .* delta ./ abs.(z))
    return out
end

function fair_l1_shrink(z, reg, delta)
    out = sign.(z) .* (abs.(z) - (delta+reg) + sqrt.((delta + reg - abs.(z)) .^2 + 4*delta .* abs(z))) ./ 2
    return out
end

function ir_gf1_shrink(z, reg, delta, a, b)
    u = a/delta
    v = b/delta
    out = sign.(z) .* (v .* abs.(z) - (1 .+ reg) + sqrt.((1 .+ reg - v.*abs.(z)).^2 + 4*(v + reg .* u) .* abs.(z))) ./ (2 * (v + reg .* u))
    return out
end

function ir_tav_shrink(z, reg, delta)
    out = zeros(size(z))
    big = reg .< abs.(z) .&& abs.(z) .< reg + delta
    out[big] = z[big] .* (1 .- reg ./ abs.(z[big]))
    big = reg + delta .<= abs.(z)
    out[big] = z[big]
    return out
end

potential_fun0 = Dict([
    (:plot, pot -> ir_potential_fun_plot(pot)),
    (:shrink, pot -> ((b, reg) -> ir_potential_fun_shrink(pot, b, reg))),
    (:potk, pot -> (z -> potk(pot,z))),
    (:wpot, pot -> (z -> wpot(pot,z))),
    (:dpot, pot -> (z -> dpot(pot,z))),
])

Base.getproperty(pot::potential, s::Symbol) =
    haskey(potential_fun0, s) ? potential_fun0[s](pot) :
    getfield(pot, s)

Base.propertynames(pot::potential) =
    (fieldnames(typeof(pot))..., keys(potential_fun0)...)

test = 0
if test == 1
    pot = potential("lange3", 10, [])
    print(pot.shrink([2,3,-5], 4))
end
