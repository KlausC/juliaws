
"""
deflation criterion.
"""
deflation_criterion{T}(sub::T, da::T, db::T) = deflation_criterion1(abs(sub), (abs(da) + abs(db) ))
@inline deflation_criterion1{T}(sub::T, da::T) = sub <= da * eps(T) / 4

