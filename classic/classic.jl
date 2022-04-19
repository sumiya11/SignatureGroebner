
BUCH_DISCARDED = 0
BUCH_NF        = 0
BUCH_SIZE      = 0
BUCH_REDUCED   = 0

function zerobuch()
    global BUCH_DISCARDED
    global BUCH_NF
    global BUCH_SIZE
    global BUCH_REDUCED
    BUCH_DISCARDED = 0
    BUCH_NF    = 0
    BUCH_SIZE      = 0
    BUCH_REDUCED   = 0
end

#-----------------------------------------------------------------------------
# Regular reductions

# try to reduce f with g in a regular way and return result
function reduction_step(f, g)
    @info "reducing $(f) with $(g)" f g

    for t in terms(f)
        # if divides
        success, u = divides(t, leading_monomial(g))
        if success
            u = divexact(t, leading_term(g))
            @info "$(leading_term(g)) | $t"

            newf = f - u*g
            @info "reduction! $(f) - $(u)*$(g) --> $newf"
            return true, newf
        end
    end
    return false, f
end

# // same as above, but in terms of set G //
# try to reduce f with G in a regular way and return result
function reduction_step(f, G::Vector{T}) where {T}
    for g in G
        success, newf = reduction_step(f, g)
        if success
            return true, newf
        end
    end
    return false, f
end

# regular normal form of f w.r.t. G
function normal_form(f, G)
    @info "computing normalform of $f w.r.t $G.."

    success = true
    newf = f
    while success
        success, newf = reduction_step(newf, G)
        @info "reduction $success"
    end
    newf
end

#-----------------------------------------------------------------------------
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

function spoly(ab)
    spoly(ab...)
end

#-----------------------------------------------------------------------------

function buchberger_criterion1(f)
    f1, f2 = f
    if isone(gcd(leading_monomial(f1), leading_monomial(f2)))
        @warn "DISCARDED $((f1, f2)) by 1st Criterion"
        global BUCH_DISCARDED
        BUCH_DISCARDED += 1
        return true
    end
    return false
end

#-----------------------------------------------------------------------------
# Groebner basis things

function groebner_basis(F)
    G = deepcopy(F)

    P = []
    for i in 1:length(G)
        for j in i+1:length(G)

            if buchberger_criterion1((G[i], G[j]))
                continue
            end

            push!(P, (G[i], G[j]))
        end
    end

    @info "generated initial G and P:"
    @info "F = $F"
    @info "G = $G"
    @info "P = $P"

    while !isempty(P)
        sort!(P, by=leading_monomial ∘ spoly)
        f = popfirst!(P)
        @info "selected $f"

        global BUCH_NF
        BUCH_NF += 1
        fNF = normal_form(spoly(f), G)

        if iszero(fNF)
            global BUCH_REDUCED
            BUCH_REDUCED += 1
        else
            for fj in G
                if buchberger_criterion1((fNF, fj)) || iszero(spoly(fNF, fj))
                    continue
                end

                push!(P, (fNF, fj))
                @info "$fNF  $fj"
                @warn "SPOLY $(last(P))"
            end

            # update G
            fNFnormalzed = fNF * inv(leading_coefficient(fNF))
            push!(G, fNFnormalzed)
        end

        @info "updated G and P"
        @info "G = $G"
        @info "P = $P"
    end

    global BUCH_SIZE
    BUCH_SIZE = length(G) - length(F)

    global BUCH_REDUCED
    global BUCH_DISCARDED

    G
end

#-----------------------------------------------------------------------------


R, (x,y,z) = PolynomialRing(QQ, ["x","y", "z"], ordering=:degrevlex)

F = [x + y + z, x*y + y*z + x*z, x*y*z - 1]

# evaluation
G = groebner_basis(F)
println("############")
println(G)
