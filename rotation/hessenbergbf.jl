"""
Module providing extensions of linalg to BigFloat and Complex{BigFloat}
"""
module LinAlgBF

import Base.LinAlg
import Base: A_mul_B!, Ac_mul_B!, A_mul_Bc!, copymutable, full
import Base.LinAlg: hessfact, hessfact!, Hessenberg, HessenbergQ
import Base.LinAlg: chkstride1, checksquare

import Base.LinAlg: reflector!, reflectorApply!
import Base.LinAlg.LAPACK: gehrd!

include("util.jl")
include("separate.jl")
include("refineprecision.jl")
include("deflationcrit.jl")
include("transformhess.jl")

typealias BigFloatOrComplex Union{Complex{BigFloat}, BigFloat}

function Hessenberg{T<:BigFloatOrComplex}(A::StridedMatrix{T})
  Hessenberg(gehrd!(A)...)
end

hessfact!{T<:BigFloatOrComplex}(A::StridedMatrix{T}) = Hessenberg(A)
hessfact{T<:BigFloatOrComplex}(A::StridedMatrix{T}) = hessfact!(copy(A))

# Perfom Householder/Hessenberg
# The Householder matrices are unitary but not hermitian in the complex case
# to provide real subdiagonal of the Hessenberg matrix.
# mimic the corresponding LAPACK.gehrd! function.
function gehrd!{T<:BigFloatOrComplex}(ilo::Integer, ihi::Integer, A::StridedMatrix{T})
  chkstride1(A)
  n = checksquare(A)
  τ = zeros(T, max(0, n - 1))
  
  ilo = max(ilo, 1)
  ihi = min(ihi, n)

  for k = ilo:ihi-1
    x = view(A, k+1:ihi, k)
    τk = reflector!(x)
    τ[k] = τk
    reflectorApply!(x, τk, view(A, k+1:n, k+1:n))
    reflectorApply!(view(A, :, k+1:n), x, τk)
  end
  A, τ
end

"""
Apply reflector from left or right, analogous to linalg/generic.jl.
"""
function reflectorApply!(A::StridedMatrix, x::AbstractVector, τ::Number)
  n, m = size(A)
  if length(x) != m
    throw(DimensionMismatch("reflector has length $(length(x)), which must match the second d    imension of matrix A, $m"))
  end
  @inbounds begin
  for j = 1:n
    # dot
    vAj = A[j,1]
    for i = 2:m
      vAj += x[i]*A[j,i]
    end

    vAj = τ*vAj

    # ger
    A[j, 1] -= vAj
    for i = 2:m
      A[j, i] -= x[i]'*vAj
    end
   end
  end
  return A
end

# various multiplications with HessenbergQ
# note the ctranspose(τ[k]) - supports the case of complex τ!
function A_mul_B!{T<:BigFloatOrComplex}(HQ::HessenbergQ{T}, A::StridedVecOrMat{T})
  n = size(A, 1)
  τ = HQ.τ
  for k = length(τ):-1:1
    τk = τ[k]'
    if τk != 0
      x = view(HQ.factors, k+1:n, k)
      reflectorApply!(x, τk, view(A, k+1:n, :))
    end
  end
  A
end

function A_mul_B!{T<:BigFloatOrComplex}(A::StridedMatrix{T}, HQ::HessenbergQ{T})
  n = size(A, 2)
  τ = HQ.τ
  for k = 1:length(τ)
    τk = τ[k]
    if τk != 0
      x = view(HQ.factors, k+1:n, k)
      reflectorApply!(view(A, :, k+1:n), x, τk)
    end
  end
  A
end

function Ac_mul_B!{T<:BigFloatOrComplex}(HQ::HessenbergQ{T}, A::StridedVecOrMat{T})
  n = size(A, 1)
  τ = HQ.τ
  for k = 1:length(τ)
    τk = τ[k]
    if τk != 0
      x = view(HQ.factors, k+1:n, k)
      reflectorApply!(x, τk, view(A, k+1:n, :))
    end
  end
  A
end

function A_mul_Bc!{T<:BigFloatOrComplex}(A::StridedMatrix{T}, HQ::HessenbergQ{T})
  n = size(A, 2)
  τ = HQ.τ
  for k = length(τ):-1:1
    τk = τ[k]'
    if τk != 0
      x = view(HQ.factors, k+1:n, k)
      reflectorApply!(view(A, :, k+1:n), x, τk)
    end
  end
  A
end

function copymutable(HQ::HessenbergQ)
  Q = eye(HQ.factors)
  A_mul_B!(HQ, Q)
  Q
end

full(HQ::HessenbergQ) = copymutable(HQ)

import Base.LinAlg.schurfact!
import Base.getindex

# Schur decomposition
immutable Schur{Ty<:BigFloatOrComplex, S<:AbstractMatrix} <: Factorization{Ty}
    T::S
    Z::S
    values::Vector
    Schur(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}, values::Vector) = new(T, Z, values)
end
Schur{Ty}(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}, values::Vector) = Schur{Ty, typeof(T)}(T, Z, values)

function getindex(F::Schur, sym::Symbol)

  if sym == :T || sym == :Schur
    F.T
  elseif sym == :Z || sym == :vectors
    F.Z
  elseif sym == :values
    F.values
  else
    error("invalid key $(sym) for Schur")
  end
end

function schurfact!{T<:BigFloatOrComplex}(A::StridedMatrix{T})
  chkstride1(A)
  n = checksquare(A)
  # transform to Hessenberg form
  HF = hessfact!(A)
  Q = full(HF[:Q])
  A = HF[:H]
  separate!(A, 1, n, Q, processPart!)
  eigv = seigendiag(A, 1, size(A, 1))
  Schur(A, Q, eigv)
end

function processPart!{T<:BigFloatOrComplex}(A::StridedMatrix{T}, ilo::Integer, ihi::Integer, Q::AbstractMatrix)
  S = ifelse(T <: Real, Float64, Complex128)
  # 1. step: estimate eigenvalues in lower precision
  F = schurfact!(S.(view(A, ilo:ihi, ilo:ihi)))
  ev = F[:values]
  refineprecision!(A, ilo, ihi, Q, ev)
end

end # module LinAlgFB

