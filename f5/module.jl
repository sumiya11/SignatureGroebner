
# import some useful things

import Base: +, -, *, <
import Base: copy, show
# provides polynomial ring and polynomials interface
using AbstractAlgebra
# provides dot product operation
import LinearAlgebra: dot

#-----------------------------------------------------------------------------

struct Signature
    index::Integer
    monom::Any
end

# our module element will store one thing: a tuple of polynomials
struct ModuleElement
    ev::Any
    sgn::Signature
end

# some magic to allow copying ModuleElement
function Base.copy(m::ModuleElement)
    return ModuleElement(deepcopy(m.ev), deepcopy(m.sgn))
end

# some magic to allow copying ModuleElement
function Base.copy(m::Signature)
    return Signature(deepcopy(m.index), deepcopy(m.monom))
end

# some magic to prin2t ModuleElement nicely
Base.show(io::IO, ::MIME"text/plain", m::ModuleElement) = println(io, "($(m.ev), $(m.sgn))")
Base.show(io::IO, m::ModuleElement) = print(io, "($(m.ev), $(m.sgn))")

#-----------------------------------------------------------------------------
# Orderings stuff

function <ₚₒₜ(a::Signature, b::Signature)
    b.index < a.index || (b.index == a.index && a.monom < b.monom)
end

>ₚₒₜ(a::Signature, b::Signature) = b <ₚₒₜ a
<(a::Signature, b::Signature) = a <ₚₒₜ b
Base.isless(a::Signature, b::Signature) = a < b

# compare two module elements w.r.t. pot strategy,
# a nice thing about julia is that we can use some latex capabilities in naming variables
function <ₚₒₜ(a::ModuleElement, b::ModuleElement)
    a.sgn <ₚₒₜ b.sgn
end

>ₚₒₜ(a::ModuleElement, b::ModuleElement) = b <ₚₒₜ a
<(a::ModuleElement, b::ModuleElement) = a <ₚₒₜ b
Base.isless(a::ModuleElement, b::ModuleElement) = a < b

#-----------------------------------------------------------------------------
# Arithmetic operations

*(t, s::Signature) = Signature(s.index, leading_monomial(s.monom*t))
*(s::Signature, t) = t * s

# returns t*m where t is something, and m is a module element
*(t, m::ModuleElement) = m * t
*(m::ModuleElement, t) = ModuleElement(m.ev*t, m.sgn*t)

# group operations
-(a::ModuleElement) = ModuleElement(-a.ev, a.sgn)
+(a::ModuleElement, b::ModuleElement) = ModuleElement(a.ev+b.ev, max(a.sgn, b.sgn))
-(a::ModuleElement, b::ModuleElement) = a + (-b)

function (≃)(a::ModuleElement, b::ModuleElement)
    a.ev == b.ev && a.sgn == b.sgn
end

#-----------------------------------------------------------------------------
# Signatures and evaluation

dot(F, m::ModuleElement) = dot(m, F)
dot(m::ModuleElement, F) = sum(m.mtuple .* F)

# evaluation hom
v(F, m::ModuleElement) = dot(F, m)
v(F, M) = [v(F, m) for m in M]

issyzygy(h::ModuleElement) = iszero(h.ev)
