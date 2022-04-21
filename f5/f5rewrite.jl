
include("LabeledPolynomial.jl")

TOTAL_DEGREES = []

debug() = false

#-----------------------------------------------------------------------------

function homogenize(F)
    base = parent(F[1])
    s = copy(symbols(base))
    push!(s, :u)
    R, xsu = PolynomialRing(base_ring(base), s, ordering=ordering(base))
    xs = xsu[1:end-1]
    u = xsu[end]

    homF = []
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

    F = []
    for homf in homF
        f = evaluate(homf, xs)
        push!(F, f)
    end

    F
end

#-----------------------------------------------------------------------------

function reduceby(f, g)
    newf = f
    for t in terms(f)
        success, u = divides(t, leading_monomial(g))
        if success
            u = divexact(t, leading_term(g))
            newf = f - u*g
            break
        end
    end
    return newf
end

function mynormalform(f, G, r)
    newf = f

    @label start
    for g in G
        newf = reduceby(newf, r[g].ev)
    end

    if newf != f
        f = newf
        @goto start
    end

    newf
end


#-----------------------------------------------------------------------------
# Groebner basis things

function f5_total_degree_lead_cmp(x, y)
    t1, t2 = total_degree(x), total_degree(y)
    if t1 < t2
        return true
    elseif t1 == t2
        return x < y
    else
        return false
    end
end

# S-polynomial

# returns u, v,  such that
# u*a - v*b is spoly(a, b)
function mults(a, b)
    u = lcm(leading_monomial(a), leading_monomial(b))
    u = divexact(u, leading_term(a))

    v = lcm(leading_monomial(a), leading_monomial(b))
    v = divexact(v, leading_term(b))

    u, v
end

# S-polynomial of a and b
function spoly(a, b)
    u, v = mults(a, b)
    u*a - v*b
end

struct RewriteRule
    index
    monom
end

function addrule!(Rule, k, sgn, r)
    s = r[k].sgn
    push!(Rule[s.index], RewriteRule(k, s.monom))
end

function findrewriting(u, k, r, Rule)
    sgn = r[k].sgn

    rules = Rule[sgn.index]
    ctr = length(rules)
    while ctr > 1
        rule = rules[ctr]
        flag, t = divides(u * sgn.monom, rule.monom)
        if flag
            @info "found rewriter of $(u * sgn) by $rule"
            return rule.index
        end
        ctr -= 1
    end

    return k
end

function isrewritable(u, k, r, Rule)
    j = findrewriting(u, k, r, Rule)
    j != k
end

function computespolys(P, r, Rule)
    S = []

    sort!(P, by= x -> x.t)

    @info "computespolys" P r Rule

    for p in P
        @info "spairs $p"
        if isrewritable(p.u, p.k, r, Rule) || isrewritable(p.v, p.l, r, Rule)
            @info "rewritten $p"
            continue
        end

        s = spoly(r[p.k].ev, r[p.l].ev)
        sgns = r[p.k].sgn * p.u
        newlp = LabeledPolynomial(s, sgns)

        push!(r, newlp)

        addrule!(Rule, length(r), sgns, r)

        @info "Added" newlp

        if !iszero(s)
            push!(S, length(r))
        end
    end

    sort!(S, by= x -> r[x].sgn)

    S
end

struct CriticalPair
    t
    k
    u
    l
    v
end


function istopreducible(f, Gprev, r)
    t = leading_monomial(f)
    for gi in Gprev
        g = r[gi].ev
        flag, _ = divides(f, leading_monomial(g))
        if flag
            @info "found top reducer" f gi r[gi]
            return true
        end
    end
    return false
end

function criticalpair(i, k, l, Gprev, r)
    tk = leading_monomial(r[k].ev)
    tl = leading_monomial(r[l].ev)

    t = lcm(tk, tl)

    u1 = divexact(t, tk)
    u2 = divexact(t, tl)

    sgn1 = r[k].sgn
    sgn2 = r[l].sgn

    @info "for lps" i k l r[k] r[l]
    @info "criticalpair reducers" sgn1 sgn2 Gprev

    if sgn1.index == i && istopreducible((u1*sgn1).monom, Gprev, r)
        @info "discarded $((k, l))" u1*sgn1
        return CriticalPair(t, 0, u1, 0, u2)
    end
    if sgn2.index == i && istopreducible((u2*sgn2).monom, Gprev, r)
        @info "discarded $((k, l))" u2*sgn2
        return CriticalPair(t, 0, u1, 0, u2)
    end

    if u1*sgn1 < u2*sgn2
        u1, u2 = u2, u1
        k, l = l, k
    end

    CriticalPair(t, k, u1, l, u2)
end

function findreductor(k, Gprev, newGcurr, r, Rule)
    t = leading_monomial(r[k].ev)
    sgnk = r[k].sgn
    for j in newGcurr
        tj = leading_monomial(r[j].ev)
        flag, u = divides(t, tj)
        if flag
            sgnj = r[j].sgn
            sgnju = sgnj*u
            if sgnju.index != sgnk.index && sgnju.monom != sgnk.monom
                if !isrewritable(u, j, r, Rule)
                    if !istopreducible(sgnju.monom, Gprev, r)
                        return j
                    end
                end
            end
        end
    end
    return 0
end

function topreduction(k, Gprev, newGcurr, r, Rule)
    if iszero(r[k].ev)
        @warn "Reduction to zero!"
        return 0, []
    end

    j = findreductor(k, Gprev, newGcurr, r, Rule)
    @info "for k=$k reductor j=$j"
    if iszero(j)
        normalize!(r[k])
        return k, []
    end

    p = r[k]
    q = r[j]
    u = divexact(leading_monomial(p.ev), leading_monomial(q.ev))
    c = leading_coefficient(p.ev) // leading_coefficient(q.ev)
    p.ev = p.ev - c*u*q.ev

    if !iszero(p.ev)
        p = normalize(p)
    end

    if u*r[j].sgn < r[k].sgn
        r[k] = p
        return 0, [k]
    else
        newr = LabeledPolynomial(p, u*r[j].sgn)
        push!(r, newr)
        addrule!(Rule, length(r), newr.sgn, r)
        return 0, [k, length(r)]
    end
end

# F4 - style reduction
function reduction(S, Gprev, Gcurr, r, Rule)
    todo = S
    completed = []

    newGcurr = copy(Gcurr)

    @info "reducing $S with"
    if debug()
        println(Gprev, " // ", Gcurr)
    end

    while !isempty(todo)
        sort!(todo, by= x -> r[x].sgn)

        k = todo[1]
        todo = todo[2:end]

        @info "reduction, handle $k, ~ $(r[k])"

        r[k].ev = mynormalform(r[k].ev, Gprev, r)
        @warn "computed normal form"
        
        newly_computed, redo = topreduction(k, Gprev, newGcurr, r, Rule)

        @info "from topreduction" newly_computed redo
        @info "todo = $todo, r = $r" r[k]

        if !iszero(newly_computed)
            push!(completed, newly_computed)
            push!(newGcurr, newly_computed)
        end

        for j in redo
            push!(todo, j)
        end

    end

    completed
end


function incremental_basis(i, Gprev, r, Rule)
    curr_idx = length(r)

    Gcurr = copy(Gprev)
    push!(Gcurr, curr_idx)

    push!(Rule, [])

    P = []
    for j in Gprev
        p = criticalpair(i, curr_idx, j, Gprev, r)
        if !iszero(p.k)
            push!(P, p)
        end
    end

    global TOTAL_DEGREES

    @info "env" Gcurr
    @info "generated pairs" P

    while !isempty(P)
        sort!(P, by= p -> total_degree(p.t))

        d = total_degree(P[1].t)
        Pd = []

        j  = 1
        while j <= length(P) && total_degree(P[j].t) == d
            push!(Pd, P[j])
            j += 1
        end
        P = P[j:end]

        @info "selected pairs" d j Pd

        S = computespolys(Pd, r, Rule)

        @info "prepared $S"
        if debug()
            println("Gcurr = ", Gcurr)
            println("r = ", r)
            println("Rule =", Rule)
        end

        R = reduction(S, Gprev, Gcurr, r, Rule)

        @info "Reduced!" R Gcurr

        for k in R
            for j in Gcurr
                p = criticalpair(i, j, k, Gprev, r)
                if p.l != 0
                    push!(P, p)
                end
            end
            push!(Gcurr, k)
        end

        if debug()
            println("###############")
            println("###############")
            println("In degree $d:")
            println("P = ", P)
            println("Gcurr = ", Gcurr)
            println("r = ", r)
            println("Rule =", Rule)
            println("###############")
            println("###############")
        end

        push!(TOTAL_DEGREES[i], d)
    end

    Gcurr
end


function signature_groebner_basis(F)
    sort!(F, by=leading_monomial, lt=f5_total_degree_lead_cmp)

    @info "after initial sort" F

    m = length(F)

    # normalize f1
    lp1 = LabeledPolynomial(F[1], Signature(1, one(F[1])))
    normalize!(lp1)

    r = [lp1]

    Gprev = [1]

    i = 2
    Rule = [[]]

    global TOTAL_DEGREES
    TOTAL_DEGREES = []
    push!(TOTAL_DEGREES, [])

    @info "at the start" Gprev Rule

    while i <= m
        push!(TOTAL_DEGREES, [])

        lpi = LabeledPolynomial(F[i], Signature(i, one(F[i])))
        normalize!(lpi)
        push!(r, lpi)

        @info "before incremental" r
        Gprev = incremental_basis(i, Gprev, r, Rule)

        if debug()
            println("Gprev =", Gprev)
        end

        i += 1
    end

    if debug()
        println("##############")
        println("##############")
    end

    [r[i].ev for i in Gprev]
end

#-----------------------------------------------------------------------------

using Logging
Logging.global_logger(ConsoleLogger(Logging.Warn))

#=
R, (x,y,z) = PolynomialRing(QQ, ["x","y", "z"], ordering=:degrevlex)

F = [x*z + 1, y*z + 1]

# evaluation
G = signature_groebner_basis(F)
Gev = [g.ev for g in G]
println("############")
println(G)
println(Gev)
=#
