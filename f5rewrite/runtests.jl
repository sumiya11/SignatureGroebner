
include("f5rewrite.jl")
using Groebner
using Test

using Logging
Logging.global_logger(Logging.ConsoleLogger(Logging.Error))

###############################################################################

rr_zerof5()

fs = Groebner.rootn(3)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))
sign_gb = [f for f in sign_gb]

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.rootn(4)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.rootn(6)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.change_ordering(Groebner.noonn(3), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.change_ordering(Groebner.noonn(4), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.change_ordering(Groebner.katsuran(4), :degrevlex)
gb = Groebner.groebner(fs)

# large
sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.change_ordering(Groebner.kinema(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
# @test Groebner.groebner(sign_gb) == gb
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.change_ordering(Groebner.sparse5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE

###############################################################################

rr_zerof5()

fs = Groebner.change_ordering(Groebner.eco5(), :degrevlex)
gb = Groebner.groebner(fs)

sign_gb = signature_groebner_basis2(homogenize(fs))

# println(sign_gb)
@test Groebner.isgroebner(sign_gb)
@error "" RR_F5_NF RR_F5_USELESS_NF RR_F5_DETECTED_USELESS_NF RR_F5_BASIS_SIZE
