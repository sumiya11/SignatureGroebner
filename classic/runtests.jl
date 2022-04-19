
include("classic.jl")
using Groebner
using Test
using AbstractAlgebra

using Logging
Logging.global_logger(Logging.ConsoleLogger(Logging.Warn))

###############################################################################

zerobuch()

fs = Groebner.rootn(3)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_NF BUCH_DISCARDED BUCH_REDUCED

###############################################################################

zerobuch()

fs = Groebner.rootn(4)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_NF BUCH_DISCARDED BUCH_REDUCED

###############################################################################

zerobuch()

fs = Groebner.change_ordering(Groebner.noonn(3), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_NF BUCH_DISCARDED BUCH_REDUCED

###############################################################################

zerobuch()

fs = Groebner.change_ordering(Groebner.noonn(4, ground=GF(2^31-1)), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_NF BUCH_DISCARDED BUCH_REDUCED

###############################################################################

zerobuch()

fs = Groebner.change_ordering(Groebner.kinema(ground=GF(2^31-1)), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_SPOLYS BUCH_DISCARDED BUCH_SPOLYS-BUCH_DISCARDED

###############################################################################

zerobuch()

fs = Groebner.change_ordering(Groebner.sparse5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_NF BUCH_DISCARDED BUCH_REDUCED

###############################################################################

zerobuch()

fs = Groebner.change_ordering(Groebner.eco5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = groebner_basis(fs)

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" BUCH_NF BUCH_DISCARDED BUCH_REDUCED
