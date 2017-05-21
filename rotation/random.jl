module RandomMatrix

export unitvector, orthogonal

using Base.LinAlg
using Base.Random

"""
return random unit vector, uniformly distributed on unit sphere of n dimensions.
"""
function unitvector{T<:AbstractFloat}(n::Integer, ::Type{T} = Float64, rng::AbstractRNG = Random.GLOBAL_RNG)

  r = zero(T)
  v = zeros(T, n)
  sqnh = sqrt(n) / 2
  while r < sqnh
    randn!(rng, v)
    r = norm(v)
  end
  v / r
end

"""
Generate random orthogonal real matrix. Determinant may be 1 or -1.
"""

function orthogonal{T<:AbstractFloat}(n::Integer, m::Integer, ::Type{T} = Float64, rng    ::AbstractRNG = Random.GLOBAL_RNG)

  Q = zeros(T,0)
  for k = 1:m
    append!(Q, unitvector(n, T, rng))
  end
  Q = reshape(Q, n, m)
  qr(Q)[1]
end

function orthogonal2{T<:AbstractFloat}(n::Integer, m::Integer, ::Type{T} = Float64, rng::AbstractRNG = Random.GLOBAL_RNG)

  n >= m || error("orthogonal random matrix requires n >= m")
  Q = eye(T, n, m)
  if m > 0 && rand(rng, Bool)
    Q[m,m] = -Q[m,m]
  end
  # R2 = zeros(T, n, m)

  for k = m:-1:1
    v = unitvector(n - k + 1, T, rng)
    # R2[k:end,k] = v
    R = UnitRotation(T)
    r = v[n-k+1]
    for i = n-k:-1:1
      g, r = givens(v[i], r, i + k - 1, i + k)
      A_mul_B!(g, R)
    end
    A_mul_B!(R', Q)
  end
  Q #, R2
end

UnitRotation{T<:AbstractFloat}(::Type{T}) = LinAlg.Rotation(LinAlg.Givens{T}[])

import Base.randn, Base.rand

function rand(rng::AbstractRNG, ::Type{BigFloat})
  bigFloat_randfill(rng, rand(rng, Float64))
end

function randn(rng::AbstractRNG, ::Type{BigFloat})
  bigFloat_randfill(rng, randn(rng, Float64))
end

"""
Append random bits in BigFloat representation of input small float number.
"""
function bigFloat_randfill{T<:Union{Float64,Float32,Float16}}(rng::AbstractRNG, a::T)
  missing = precision(BigFloat) - precision(a)
  bigi = BigInt(1) << missing
  add = ( rand(rng, BigInt(0):bigi-1) / bigi )
  add *= big"2.0"^(exponent(a) - precision(a))
  add = copysign(add, a)
  a + add
end

end


