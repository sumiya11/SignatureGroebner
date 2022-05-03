
include("LabeledPolynomial.jl")

TOTAL_DEGREES = []

TOPREDUCTIONSF5 = 0
NORMALFORMS = 0

RR_F5_DETECTED_USELESS_NF = 0
RR_F5_USELESS_NF   = 0
RR_F5_NF        = 0
RR_F5_BASIS_SIZE      = 0

function rr_zerof5()
    global RR_F5_DETECTED_USELESS_NF
    global RR_F5_USELESS_NF
    global RR_F5_NF
    global RR_F5_BASIS_SIZE
    RR_F5_DETECTED_USELESS_NF = 0
    RR_F5_USELESS_NF        = 0
    RR_F5_NF      = 0
    RR_F5_BASIS_SIZE   = 0
end

function katsuran(n; ground=QQ)
    R, x = PolynomialRing(ground, ["x$i" for i in 0:n])

    return [
        (sum(x[abs(l)+1]*x[abs(m-l)+1] for l=-n:n if abs(m-l)<=n) -
        x[m+1] for m=0:n-1)...,
        x[1] + 2sum(x[i+1] for i=1:n) - 1
    ]
end

using OffsetArrays

debug() = false

#-----------------------------------------------------------------------------

without_k(arr, k) = arr[1:length(arr) .!= k]

function interreduce_basis(F::Vector{AbstractAlgebra.Generic.MPoly{T}}) where {T}

    F = Vector{AbstractAlgebra.Generic.MPoly{T}}(F)

    @label start

    newF = empty(F)

    for (i, f) in enumerate(F)
        reducers = union(newF, F[i+1:end])
        _, nf = divrem(f, reducers)
        if !iszero(nf)
            push!(newF, nf)
        end
    end

    if F != newF
        F = newF
        @goto start
    end

    newF
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

function reduceby(f, g)
    newf = f
    leadg = leading_monomial(g)
    leadf = leading_monomial(f)
    success, u = divides(leadf, leadg)
    if success
        u = divexact(leadf, leadg)
        newf = f - u*g
    end
    return success, newf
end

function mynormalform(f, G, r)

    global NORMALFORMS
    NORMALFORMS += 1

    a, b = divrem(f, [r[g].ev for g in G])
    return b

    newf = f

    @label start
    updated = false
    for g in G
        success, newf = reduceby(newf, r[g].ev)
        updated = updated || success
    end

    if updated
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

struct RewriteRule{Poly}
    index::Int
    monom::Poly
end

function addrule!(Rule::Vector{Vector{RewriteRule{Poly}}}, sgn, k) where {Poly}
    push!(Rule[sgn.index], RewriteRule{Poly}(k, sgn.monom))
end

function findrewriting(u, k, r, Rule)
    sgn = r[k].sgn

    rules = Rule[sgn.index]
    ctr = length(rules)

    debug() && @info "search for rewriter $rules" u k sgn*u

    while ctr >= 1
        rule = rules[ctr]
        flag, t = divides(u * sgn.monom, rule.monom)
        if flag
            debug() && @info "found rewriter of $(u * sgn) by $rule"
            return rule.index
        end
        ctr -= 1
    end

    return k
end

# changed
function findrewriting2(u, k, r, Rule)
    sgn = r[k].sgn

    rules = Rule[sgn.index]
    ctr = length(rules)

    @info "search for rewriter $rules" u k sgn*u

    println(rules)
    println("EEEE")

    while ctr >= 1
        rule = rules[ctr]
        flag, t = divides(u * sgn.monom, rule.monom)
        if flag
            debug() && @info "found rewriter of $(u * sgn) by $rule"
            if rule.index != k
                return rule.index
            end
        end
        ctr -= 1
    end

    return k
end

function isrewritable(u, k, r, Rule)
    j = findrewriting(u, k, r, Rule)
    j != k
end

function computespolys(P, r::Vector{Poly}, Rule) where {Poly}
    S = Int[]

    sort!(P, by= x -> x.t)

    debug() && @info "computespolys" P r Rule

    for p in P
        debug() && @info "spairs $p"

        # ADDED
        sgnk = r[p.k].sgn
        sgnl = r[p.l].sgn
        debug() && @info "is rewritable?? $(p.u) or $(p.v)"
        if isrewritable(p.u, p.k, r, Rule) || isrewritable(p.v, p.l, r, Rule)
            global RR_F5_DETECTED_USELESS_NF
            RR_F5_DETECTED_USELESS_NF += 1
            
            debug() && @info "rewritten $p"
            continue
        end

        s = spoly(r[p.k].ev, r[p.l].ev)
        sgns = r[p.k].sgn * p.u
        newlp = LabeledPolynomial(s, sgns)

        push!(r, newlp)

        addrule!(Rule, sgns, length(r))

        # @warn "Added" newlp
        # println(r)

        if !iszero(s)
            push!(S, length(r))
        end
    end

    sort!(S, by= x -> r[x].sgn)

    S
end

struct CriticalPair{Poly}
    t::Poly
    k::Int
    u::Poly
    l::Int
    v::Poly
end

# some magic to prin2t ModuleElement nicely
Base.show(io::IO, ::MIME"text/plain", m::CriticalPair) = println(io, "$((m.t, m.k, m.u, m.l, m.v))")
Base.show(io::IO, m::CriticalPair) = print(io, "$((m.t, m.k, m.u, m.l, m.v))")


function istopreducible(f, Gprev, r)
    t = leading_monomial(f)
    debug() && @info "istopreducible $f by " Gprev r[Gprev]
    for gi in Gprev
        g = r[gi].ev
        flag, _ = divides(t, leading_monomial(g))
        if flag
            debug() && @info "found top reducer" f gi r[gi]
            return true
        end
    end
    return false
end

function criticalpair(i, k, l, Gprev, r, Rule)
    tk = leading_monomial(r[k].ev)
    tl = leading_monomial(r[l].ev)

    t = lcm(tk, tl)

    u1 = divexact(t, tk)
    u2 = divexact(t, tl)

    sgn1 = r[k].sgn
    sgn2 = r[l].sgn

    # debug() && @info "for lps" i k l r[k] r[l]
    # debug() && @info "criticalpair reducers" sgn1 sgn2 Gprev

    debug() && @info "for lps" i k l r[k] r[l]
    debug() && @info "criticalpair reducers" sgn1 sgn2 Gprev

    usgn1 = u1*sgn1
    usgn2 = u2*sgn2

    if usgn1 == usgn2
        global RR_F5_DETECTED_USELESS_NF
        RR_F5_DETECTED_USELESS_NF += 1

        return CriticalPair(t, 0, u1, 0, u2)
    end

    @assert usgn2 != usgn1

    debug() && @info "multiplied sgns" usgn1 usgn2

    if sgn1.index == i && istopreducible(usgn1.monom, Gprev, r)
        global RR_F5_DETECTED_USELESS_NF
        RR_F5_DETECTED_USELESS_NF += 1

        debug() && @info "discarded $((k, l))" u1*sgn1
        return CriticalPair(t, 0, u1, 0, u2)
    end
    if sgn2.index == i && istopreducible(usgn2.monom, Gprev, r)
        global RR_F5_DETECTED_USELESS_NF
        RR_F5_DETECTED_USELESS_NF += 1

        debug() && @info "discarded $((k, l))" u2*sgn2
        return CriticalPair(t, 0, u1, 0, u2)
    end

    # @error "is??" usgn1.monom, l, r[l]

    # Added
    # if istopreducible(usgn1.monom, [l], r)
    #    @error "NEW discarded $((k, l))" u2*sgn2
    #    return CriticalPair(t, 0, u1, 0, u2)
    # end

    if isrewritable(u1, k, r, Rule) || isrewritable(u2, l, r, Rule)
        global RR_F5_DETECTED_USELESS_NF
        RR_F5_DETECTED_USELESS_NF += 1

        debug() && @info "rewritten $u1 or $u2"
        return CriticalPair(t, 0, u1, 0, u2)
    end

    if u1*sgn1 < u2*sgn2
        u1, u2 = u2, u1
        k, l = l, k
    end

    # debug() && @info "CriticalPair" t k u1 l u2
    debug() && @info "CriticalPair" t k u1 l u2

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
            if sgnju.index != sgnk.index || sgnju.monom != sgnk.monom
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
        @error "Reduction to zero!"

        global RR_F5_USELESS_NF
        RR_F5_USELESS_NF += 1

        return 0, ()
    end

    j = findreductor(k, Gprev, newGcurr, r, Rule)

    (debug() || false) && @info "for k=$k reductor j=$j"

    if iszero(j)
        normalize!(r[k])
        return k, ()
    end

    p = r[k]
    q = r[j]

    (debug() || false) && println("reductor poly $q")

    u = divexact(leading_monomial(p.ev), leading_monomial(q.ev))
    c = leading_coefficient(p.ev) // leading_coefficient(q.ev)

    global TOPREDUCTIONSF5
    TOPREDUCTIONSF5 += 1

    newpev = p.ev - c*u*q.ev

    (debug() || false) && println("reductor res $p")

    if !iszero(newpev)
        newpev = normalize(newpev)
    end

    if u*r[j].sgn < r[k].sgn
        r[k].ev = newpev
        return 0, (k,)
    else
        newr = LabeledPolynomial(newpev, u*r[j].sgn)
        push!(r, newr)
        addrule!(Rule, newr.sgn, length(r))
        return 0, (k, length(r))
    end
end

function reduction(S, Gprev, Gcurr, r, Rule)
    todo = S
    completed = Int[]

    newGcurr = copy(Gcurr)

    debug() && @info "reducing $S with"
    if debug()
        println(Gprev, " // ", Gcurr)
    end

    if debug() || false
        println("r = $r")
    end

    while !isempty(todo)
        sort!(todo, by= x -> r[x].sgn)

        k = todo[1]
        todo = todo[2:end]

        (debug() || false) && @info "reduction, handle $k ($(r[k])) by $(r[Gprev])"

        global RR_F5_NF
        RR_F5_NF += 1
        r[k].ev = mynormalform(r[k].ev, Gprev, r)

        (debug() || false) && @info "from normal form" r[k]

        if debug() || false
            println("GOING TO TOPREDUCTION with Gprev = $Gprev and Gcurr=$newGcurr")
        end

        newly_computed, redo = topreduction(k, Gprev, newGcurr, r, Rule)

        (debug() || false) && @info "After top" r[k]

        (debug() || false) && @info "from topreduction" newly_computed redo
        (debug() || false) && @info "todo = $todo"

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


function incremental_basis(i, Gprev, r::Vector{LabeledPolynomial{Poly}}, Rule) where {Poly}
    curr_idx = length(r)

    Gcurr = copy(Gprev)
    push!(Gcurr, curr_idx)

    push!(Rule, Vector{RewriteRule{Poly}}(undef, 0) )

    P = Vector{CriticalPair{Poly}}(undef, 0)
    for j in Gprev
        p = criticalpair(i, curr_idx, j, Gprev, r, Rule)
        if !iszero(p.k)
            push!(P, p)
        end
    end

    # global TOTAL_DEGREES

    debug() && @info "env" Gcurr
    @warn "$i iter, generated pairs" length(P) # P

    if debug() || false
        println("Pairs are:")
        for p in P
            println("($(p.k), $(p.l))")
        end
    end

    while !isempty(P)
        sort!(P, by= p -> total_degree(p.t))

        d = total_degree(P[1].t)
        Pd = Vector{CriticalPair{Poly}}(undef, 0)

        j  = 1
        while j <= length(P) && total_degree(P[j].t) == d
            push!(Pd, P[j])
            j += 1
        end
        # j = max(1, div(j, 5))
        # Pd = Pd[1:j]
        P = P[j:end]

        @warn "SELECTED $(length(Pd)) pairs of degree $d out of $(length(Pd) + length(P))"

        # println("Spolys: $Pd")
        # println("$(Pd[1].k)-th poly $(r[Pd[1].k])")
        # println("$(Pd[1].l)-th poly $(r[Pd[1].l])")

        if length(Pd) > 1
            # println("$(Pd[2].k)-th poly $(r[Pd[2].k])")
            # println("$(Pd[2].l)-th poly $(r[Pd[2].l])")
        end

        debug() && @info "selected pairs" d j Pd

        S = computespolys(Pd, r, Rule)

        if (debug() || false) && (length(S) > 0)
            println("# SPOLYS $(length(S))")
            println(S)
            println("First spoly: $(r[S[1]])")
        end

        R = reduction(S, Gprev, Gcurr, r, Rule)

        if (debug() || false) && (length(R) > 0)
            println("AFTER reduction $(length(R))")
            println(R)
            println("First reduced: $(r[R[1]])")
        end

        for k in R
            for j in Gcurr
                p = criticalpair(i, j, k, Gprev, r, Rule)
                if p.l != 0
                    push!(P, p)
                end
            end
            (debug() || true) && println("New # pairs for poly $k:", length(P))
            push!(Gcurr, k)
            (debug() || true) &&  println("New # Gcurr:", length(Gcurr))
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

        # push!(TOTAL_DEGREES[i], d)
    end

    if debug() || false
        println("ENDED with $Gcurr polynomials")
    end

    Gcurr
end

function getreducers(i, Gprev, r)
    ans = Gprev[1:(i-1)]
    ans = filter!(x -> !iszero(r[x].ev), append!(ans, Gprev[(i+1):end]))
    ans
end

function interreduction(Gprev, r)

    @label start
    updated = false

    for (i, idx) in enumerate(Gprev)
        reducers = getreducers(i, Gprev, r)
        newf = mynormalform(r[idx].ev, reducers, r)
        # @warn "" i idx newf r[idx].ev
        updated = updated || (newf != r[idx].ev)
        r[idx].ev = newf
    end

    if updated
        @goto start
    end

    res = []
    for i in Gprev
        if !iszero(r[i].ev)
            push!(res, r[i].ev)
        end
    end

    return res
end

function setup_reduced_basis(Gprev, r, Rule)
    B = interreduction(Gprev, r)

    Gcurr = collect(1:length(B))

    resize!(r, length(B))
    r[1:length(B)] = [
        LabeledPolynomial(B[i], Signature2(i, one(B[i])))
        for i in 1:length(B)
    ]

    resize!(Rule, length(B)+1)
    for i in 0:length(Rule)-1
        Rule[i] = []
    end

    for j in 1:length(B)-1
        t = leading_monomial(B[j])
        for k in j + 1:length(B)
            u = divexact(lcm(t, leading_monomial(B[k])), leading_monomial(B[k]))
            addrule!(Rule, Signature2(k, u), 0)
        end
    end

    return Gcurr
end

function signature_groebner_basis2(F::Vector{Poly}) where {Poly}
    sort!(F, by=leading_monomial, lt=f5_total_degree_lead_cmp)

    debug() && @warn "after initial sort" F

    # F = interreduce_basis(F)

    debug() && @warn "after interreduction" F

    m = length(F)

    # normalize f1
    lp1 = LabeledPolynomial(F[1], Signature2(1, one(F[1])))
    normalize!(lp1)

    r = [lp1]

    Gprev = [1]

    i = 2
    Rule = Vector{Vector{RewriteRule{Poly}}}(undef, 0)
    push!(Rule, Vector{RewriteRule{Poly}}(undef, 0))
    push!(Rule, Vector{RewriteRule{Poly}}(undef, 0))
    # Rule = OffsetArray(Rule, 0:1)

    # global TOTAL_DEGREES
    # TOTAL_DEGREES = []
    # push!(TOTAL_DEGREES, [])

    @debug "at the start" Gprev Rule

    while i <= m
        # push!(TOTAL_DEGREES, [])

        lpi = LabeledPolynomial(F[i], Signature2(i, one(F[i])))
        normalize!(lpi)
        push!(r, lpi)


        @debug "before incremental" r
        Gprev = incremental_basis(i, Gprev, r, Rule)

        # Gprev = setup_reduced_basis(Gprev, r, Rule)

        #println(Gprev)
        #println(r)
        #println(Rule)

        if debug()
            println("Gprev =", Gprev, r, Rule)
        end

        i += 1
    end

    if debug()
        println("##############")
        println("##############")
    end

    global RR_F5_BASIS_SIZE
    RR_F5_BASIS_SIZE = length(Gprev)

    # println(Gprev)
    # println(r[Gprev])
    B = interreduction(Gprev, r)

    sort(B, by=leading_monomial, rev=true)
end

#-----------------------------------------------------------------------------

using Logging
Logging.global_logger(ConsoleLogger(Logging.Info))

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
