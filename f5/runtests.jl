
include("f5.jl")
using Groebner
using Test

using Logging
Logging.global_logger(Logging.ConsoleLogger(Logging.Error))

###############################################################################

zerof5()

fs = Groebner.rootn(3)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.rootn(4)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.rootn(6)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.noonn(3), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.noonn(4), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.katsura6(), :degrevlex)
gb = Groebner.groebner(fs)

# large
sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.kinema(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.sparse5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.eco5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(fs)
sign_gb = [f.ev for f in sign_gb]

println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_DISCARDED F5_REDUCED
