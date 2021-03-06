
include("module.jl")

using Logging
Logging.global_logger(Logging.ConsoleLogger(Logging.Info))

F5_DETECTED_USELESS_NF = 0
F5_USELESS_NF   = 0
F5_NF        = 0
F5_BASIS_SIZE      = 0

function zerof5()
    global F5_DETECTED_USELESS_NF
    global F5_USELESS_NF
    global F5_NF
    global F5_BASIS_SIZE
    F5_DETECTED_USELESS_NF = 0
    F5_USELESS_NF        = 0
    F5_NF      = 0
    F5_BASIS_SIZE   = 0
end


function katsuran(n; ground=QQ)
    R, x = PolynomialRing(ground, ["x$i" for i in 0:n])

    return [
        (sum(x[abs(l)+1]*x[abs(m-l)+1] for l=-n:n if abs(m-l)<=n) -
        x[m+1] for m=0:n-1)...,
        x[1] + 2sum(x[i+1] for i=1:n) - 1
    ]
end


#-----------------------------------------------------------------------------

function homogenize(F)
    base = parent(F[1])
    s = copy(symbols(base))
    push!(s, :h)
    R, xsu = PolynomialRing(base_ring(base), s, ordering=ordering(base))
    xs = xsu[1:end-1]
    u = xsu[end]

    homF = Vector{eltype(F)}(undef, 0)
    for f in F
        tdeg = total_degree(leading_monomial(f))
        f = evaluate(f, xs)
        homf = zero(R)
        for t in terms(f)
            homf += t * u^(tdeg - total_degree(t))
        end
        push!(homF, homf)
    end

    homF
end

function dehomogenize(homF)
    R = parent(homF[1])
    xs = collect(gens(R))
    xs[end] = one(R)

    F = Vector{eltype(homF)}(undef, 0)
    for homf in homF
        f = evaluate(homf, xs)
        push!(F, f)
    end

    F
end

#-----------------------------------------------------------------------------
# Regular reductions

# try to reduce f with g in a regular way and return result
function regular_reduction_step(f::ModuleElement, g::ModuleElement)
    evalf = f.ev
    evalg = g.ev

    # @info "reducing $(f) with $(g)" evalf evalg

    for t in terms(evalf)
        # if divides
        success, u = divides(t, leading_term(evalg))
        if success
            u = divexact(t, leading_term(evalg))
            @info "$(leading_term(evalg)) | $t"

            # if reduction is regular (and does not change signature)
            # not really correct!
            if f > u*g
                newf = f - u*g
                # @info "regular! $(f) - $(u)*$(g) --> $newf"
                return true, newf
            end
        end
    end
    return false, f
end

# // same as above, but in terms of set G //
# try to reduce f with G in a regular way and return result
function regular_reduction_step(f::ModuleElement, G)
    evalf = f.ev
    for g in G
        success, newf = regular_reduction_step(f, g)
        if success
            return true, newf
        end
    end
    return false, f
end

# regular normal form of f w.r.t. G
function regular_normal_form(f::ModuleElement, G)
    # @info "computing reg normalform of $f w.r.t $G.."

    success = true
    newf = copy(f)
    while success
        success, newf = regular_reduction_step(newf, G)
        # @info "reduction $success"
    end
    newf
end

#-----------------------------------------------------------------------------
# Singular reductions

# if f is singurarly top reducible by G
function issingularlytopreducible(f::ModuleElement, G)
    leadevalf = leading_monomial(f.ev)
    for g in G
        evalg = g.ev
        success, u = divides(leadevalf, leading_monomial(evalg))
        if success
            if (g*u).sgn == f.sgn
                return true
            end
        end
    end
    return false
end

# if f is singurarly top reducible by G
function syzygy_criterion(f::ModuleElement, syzygies)
    for syz in syzygies
        success, u = divides(f.sgn.monom, syz.sgn.monom)
        if success && f.sgn.index == syz.sgn.index
            # @warn "DISCARDED $f by $syz"
            global F5_DETECTED_USELESS_NF
            F5_DETECTED_USELESS_NF += 1
            return true
        end
    end
    return false
end

#-----------------------------------------------------------------------------
# S-polynomial

# returns u, v,  such that
# u*a - v*b is spoly(a, b)
function mults(a::ModuleElement, b::ModuleElement)
    u = lcm(leading_monomial(a.ev), leading_monomial(b.ev))
    u = divexact(u, leading_term(a.ev))

    v = lcm(leading_monomial(a.ev), leading_monomial(b.ev))
    v = divexact(v, leading_term(b.ev))

    u, v
end

# S-polynomial of a and b
function spoly(a::ModuleElement, b::ModuleElement)
    u, v = mults(a, b)
    u*a - v*b
end

#-----------------------------------------------------------------------------
# Groebner basis things

function construct_module(F)
    [
        ModuleElement(f, Signature(i, one(f)))
        for (i, f) in enumerate(F)
    ]
end

function signature_groebner_basis(F)
    sort!(F, by=leading_monomial, lt=f5_total_degree_lead_cmp)

    F = map(f -> map_coefficients(c -> c // leading_coefficient(f), f), F)

    G = construct_module(F)

    syzygies = []
    for i in 1:length(G)
        for j in i+1:length(G)
            s1 = (G[i]*leading_monomial(G[j].ev)).sgn
            s2 = (G[j]*leading_monomial(G[i].ev)).sgn
            push!(syzygies, ModuleElement( zero(F[1]), max(s1, s2) ) )
        end
    end

    P = []
    for i in 1:length(G)
        for j in i+1:length(G)
            if syzygy_criterion(spoly(G[i], G[j]), syzygies)
                continue
            end

            push!(P, spoly(G[i], G[j]))
        end
    end

    # @info "generated initial G and P:"
    # @info "F = $F"
    # @info "G = $G"
    # @info "P = $P"
    # @info "syz = $syzygies"

    while !isempty(P)
        sort!(P, rev=false)
        f = popfirst!(P)
        @info "selected $f"

        if syzygy_criterion(f, syzygies)
            continue
        end

        # global F5_NF
        # F5_NF += 1
        fNF = regular_normal_form(f, G)
        global F5_NF
        F5_NF += 1

        if issyzygy(fNF)
            # @warn "Reduction to zero!"
            # global F5_REDUCED
            # F5_REDUCED += 1
            push!(syzygies, fNF)

            global F5_USELESS_NF
            F5_USELESS_NF += 1
        elseif !issingularlytopreducible(fNF, G)
            # update P
            for fj in G

                if syzygy_criterion(spoly(fNF, fj), syzygies)
                    continue
                end

                push!(P, spoly(fNF, fj))

                # @info "$fNF  $fj"
                # @info "SPOLY $(last(P))"
            end

            # update G
            fNFnormalzed = fNF * inv(leading_coefficient(fNF.ev))
            push!(G, fNFnormalzed)
        else
            global F5_USELESS_NF
            F5_USELESS_NF += 1
        end

        # @info "updated G and P"
        # @info "G = $G"
        # @info "P = $P"
        # @info "syz = $syzygies"
    end

    # global F5_SIZE
    # F5_SIZE = length(G) - length(F)

    # global F5_REDUCED
    # global F5_DISCARDED
    global F5_BASIS_SIZE
    F5_BASIS_SIZE = length(G)

    G
end

#-----------------------------------------------------------------------------

using Logging
Logging.global_logger(ConsoleLogger(Logging.Warn))

R, (x,y,z) = PolynomialRing(QQ, ["x","y", "z"], ordering=:degrevlex)

F = [x*z + 1, y*z + 1]

# evaluation
G = signature_groebner_basis(F)
Gev = [g.ev for g in G]
println("############")
println(G)
println(Gev)
