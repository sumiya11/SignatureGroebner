
# import some useful things

import Base: +, -, *, <, ==
import Base: copy, show
# provides polynomial ring and polynomials interface
using AbstractAlgebra
# provides dot product operation
import LinearAlgebra: dot

#-----------------------------------------------------------------------------
# Signature2 stuff

struct Signature2{Poly}
    index::Int
    monom::Poly
end

function <ₛ(s1::Signature2, s2::Signature2)
    s1.index < s2.index || (s1.index == s2.index && s1.monom < s2.monom)
end

==(a::Signature2, b::Signature2) = a.index == b.index && a.monom == b.monom

<(a::Signature2, b::Signature2) = a <ₛ b
Base.isless(a::Signature2, b::Signature2) = a < b

#-----------------------------------------------------------------------------

# our module element will store one thing: a tuple of polynomials
mutable struct LabeledPolynomial{Poly}
    ev::Poly
    sgn::Signature2{Poly}
end

function normalize(f::LabeledPolynomial)
    lead = leading_coefficient(f.ev)
    LabeledPolynomial(map_coefficients(c -> c // lead, f.ev), f.sgn)
end

function normalize(f)
    lead = leading_coefficient(f)
    map_coefficients(c -> c // lead, f)
end

function normalize!(f::LabeledPolynomial)
    lead = leading_coefficient(f.ev)
    f.ev = map_coefficients(c -> c // lead, f.ev)
    f
end

# some magic to prin2t ModuleElement nicely
Base.show(io::IO, ::MIME"text/plain", m::LabeledPolynomial) = println(io, "($(m.ev), $(m.sgn))")
Base.show(io::IO, m::LabeledPolynomial) = print(io, "($(m.ev), $(m.sgn))")

Base.show(io::IO, m::Signature2) = print(io, "<$(m.index), $(m.monom)>")
Base.show(io::IO, ::MIME"text/plain", m::Signature2) = println(io, "<$(m.index), $(m.monom)>")

#-----------------------------------------------------------------------------
# Arithmetic operations

*(t::Poly, s::Signature2{Poly}) where {Poly} = Signature2(s.index, s.monom*t)
*(s::Signature2{Poly}, t::Poly) where {Poly} = t * s
