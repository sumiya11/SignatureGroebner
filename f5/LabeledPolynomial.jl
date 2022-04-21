
# import some useful things

import Base: +, -, *, <
import Base: copy, show
# provides polynomial ring and polynomials interface
using AbstractAlgebra
# provides dot product operation
import LinearAlgebra: dot

#-----------------------------------------------------------------------------
# Signature stuff

struct Signature{Poly}
    index::Int
    monom::Poly
end

function <ₛ(s1::Signature, s2::Signature)
    s1.index < s2.index || (s1.index == s2.index && s1.monom < s2.monom)
end

<(a::Signature, b::Signature) = a <ₛ b
Base.isless(a::Signature, b::Signature) = a < b

#-----------------------------------------------------------------------------

# our module element will store one thing: a tuple of polynomials
mutable struct LabeledPolynomial{Poly}
    ev::Poly
    sgn::Signature{Poly}
end

function normalize(f::LabeledPolynomial)
    lead = leading_coefficient(f.ev)
    LabeledPolynomial(map_coefficients(c -> c // lead, f.ev), f.sgn)
end

function normalize!(f::LabeledPolynomial)
    lead = leading_coefficient(f.ev)
    f.ev = map_coefficients(c -> c // lead, f.ev)
    f
end

# some magic to prin2t ModuleElement nicely
Base.show(io::IO, ::MIME"text/plain", m::LabeledPolynomial) = println(io, "($(m.ev), $(m.sgn))")
Base.show(io::IO, m::LabeledPolynomial) = print(io, "($(m.ev), $(m.sgn))")

Base.show(io::IO, m::Signature) = print(io, "<$(m.index), $(m.monom)>")
Base.show(io::IO, ::MIME"text/plain", m::Signature) = println(io, "<$(m.index), $(m.monom)>")

#-----------------------------------------------------------------------------
# Arithmetic operations

*(t::Poly, s::Signature{Poly}) where {Poly} = Signature(s.index, s.monom*t)
*(s::Signature{Poly}, t::Poly) where {Poly} = t * s
