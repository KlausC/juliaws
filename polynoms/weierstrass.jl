"""
Approximate x -> abs(x) by ordinary polynomials on [-1,1]
"""
module Weierstrass

using Polynomials

function absapprox{T<:Real}(::Type{T}, n::Integer)

  # Polynom T[1] + x * ( T[2] + x * ... )
  x = Poly(T[0, 0, 1])
  a = Poly(T[0])
  for k = 1:n
    a = (x - a * a)  + a + a
    a.a .//= 2
  end
  a
end

function absapprox{T<:Real}(x::T, n::Integer)

  # Polynom T[1] + x * ( T[2] + x * ... )
  x = x ^ 2
  a = zero(T)
  for k = 1:n
    a = (x - a * a) // 2  + a
  end
  a
end

function absapprox{T<:Real}(x::AbstractVector{T}, n::Integer)

  # Polynom T[1] + x * ( T[2] + x * ... )
  x = x .^ 2
  a = zeros(x)
  for k = 1:n
    a .= (x - a .* a) / 2  + a
  end
  a
end

end # module

