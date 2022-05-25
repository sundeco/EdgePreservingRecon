export potential

mutable struct potential
    ptype::String
    delta
    param
end

function potential_fun(ptype::String, delta ; param = [])
    scale = 1
    return potential(ptype, delta, param)
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
    else
        error("not done yet")
    end
end

function wpot(pot::potential, z)
    if pot.ptype == "lange3"
        return 1 ./ (1 .+ abs.(z ./ pot.delta))
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
    else
        error("not done yet")
    end
end

function dpot(pot::potential, z)
    return z ./ (1 .+ abs.(z ./ pot.delta))
end

function ir_potential_fun_plot(pot::potential)
    error("not done yet")
    z = LinRange(-1, 1, 101) * 2
    z = z * pot.delta
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
