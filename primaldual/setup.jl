# julia
#
include("problem.jl")
#include("pdutils.jl")
using MPS
using Base.SparseArrays.SPQR
#ios = openmps("x")
#p1 = LPData()
#readmps!(p1, ios)

#p2 = LPData()
#ios2 = openmps("data/maros-r7")
#readmps!(p2, ios2)

function problem(name::AbstractString)
  p = LPData()
  ios = openmps("data/" * name)
  readmps!(p, ios)
  p
end
