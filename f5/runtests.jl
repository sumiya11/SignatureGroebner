
include("f5.jl")
using Groebner
using Test

using Logging
Logging.global_logger(Logging.ConsoleLogger(Logging.Error))

###############################################################################

zerof5()

fs = Groebner.rootn(3)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.rootn(4)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.rootn(6)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.noonn(3), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.noonn(4), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.change_ordering(katsuran(4), :degrevlex)
gb = Groebner.groebner(fs)

# large
sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.kinema(ground=GF(2^31-1)), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.sparse5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE

###############################################################################

zerof5()

fs = Groebner.change_ordering(Groebner.eco5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis(homogenize(fs))
sign_gb = [f.ev for f in sign_gb]

# println(sign_gb)
# @test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" F5_NF F5_USELESS_NF F5_DETECTED_USELESS_NF F5_BASIS_SIZE
