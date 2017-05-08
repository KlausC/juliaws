module pythagorean


# construct random coprime pythagorean triple a^2 + b^2 = c^2
# 
#
function rand_coprimes(bits::Int)
  maxi = BigInt(2) ^ (bits-2) * 2 + 1
  req = 2 * maxi * (maxi - 1)
  T = mintype(req)
  r = 1:T(req)
  p = T(1) // T(1)
  while isodd(p.den) && isodd(p.num)
    p = rand(r) // rand(r)
  end
  extrema((p.den, p.num))
end


# uniform distribution with respect to angles in [0,pi/2] by atan(b / a)
# delivers coprime m and n of different parity.
# The coarsity of the distribution is controlled by the real tol-value >= 0.0
# used in rationalize.
#
function rand_coprimes(tol::Float64)
  f(s::Float64) = s / ( 1.0 + sqrt(1 - s * s) )
  p = 1 // 1
  while p.den == 1 # isodd(p.den) && isodd(p.num)
    p = rationalize(BigInt, f(sin(rand()*pi/2)), tol=tol)
  end
  extrema((p.den, p.num))
end

# determine the smallest signed integer type, which is capable of representing the given value
function mintype{S<:Integer,T<:Integer,U<:Integer}(value::S, fact1::T = 1, fact2::U = 1)
  value == 0 && return Int8
  tmaxmin = value > 0 ? typemax : typemin 
  for TT in (Int8, Int16, Int32, Int64, Int128)
    TT(1) <= tmaxmin(TT) ÷ value ÷ fact1 ÷ fact2 && return TT
  end
  BigInt
end

# determine the smallest unsigned integer type, which is capable of representing the given value
function mintype{S<:Unsigned,T<:Unsigned,U<:Unsigned}(value::S, fact1::T = 0x1, fact2::U = 0x1)
  value == 0 && return UInt8
  for TT in (UInt8, UInt16, UInt32, UInt64, UInt128)
    TT(1) <= typemax(TT) ÷ value ÷ fact1 ÷ fact2 && return TT
  end
  BigInt
end

#convert integer to integer of minimal type with same value.
mintyped(n::Integer) = mintype(n)(n)

function rand_pythagorean_triples(bits)
  tmnabc(rand_coprimes(bits)...)
end

# transform m,n - form to a,b,c - form
function tmnabc{T<:Integer, S<:Integer}(m::T, n::S)
  T1 = promote_type(T, S)
  m, n = T1(m), T1(n)
  TT = mintype(m, m, T1(2))
  m2, n2 = TT(m)^2, TT(n)^2
  m2 - n2, TT(m) * TT(n) * TT(2), m2 + n2
end

# transform a,b,c - form to m,n - formrgu
function tabcmn{T<:Integer, S<:Integer, U<:Integer}(a::T, b::S, c::U)
  T1 = promote_type(T, S, U)
  cca = gcd(c, a)
  ( gcd(c, b) == cca ) || throw(ArgumentError("no common gcd"))
  c ÷= cca
  a ÷= cca
  b ÷= cca
  c >= a || begin a, c = c, a; end
  c >= b || begin a, c = c, b; end
  isodd(c) || throw(ArgumentError("c must be odd"))
  ( isodd(a) != isodd(b) ) || throw(ArgumentError("exactly one of a, b must be odd"))
  if isodd(b)
    a, b = b, a
  end
  m2 = (c + a) ÷ 2
  n2 = (c - a) ÷ 2
  m, n = isqrt(m2), isqrt(n2)
  ( m*m == m2 ) || throw(ArgumentError("half sum of odds is not square"))
  ( n*n == n2 ) || throw(ArgumentError("half diff of odds is not square"))
  ( m * n * 2 == b ) || throw(ArgumentError("sum of squares of a and is not square of c"))
  m, n, cca
end

tabcmn(t::Tuple) = length(t) == 3 ? tabcmn(t...) : throw(ArgumentError("tuple size 3 required"))
tmnabc(t::Tuple) = length(t) == 2 ? tmnabc(t...) : throw(ArgumentError("tuple size 2 required"))

# transform generating Matrix form 2x2 to 3x3 form
function tmnabc_genmat{T<:Integer}(U::AbstractArray{T,2})
  size(U) = (2, 2) || throw(ArgumentError("2x2 matrix required"))
  A = similar(U, (3, 3))
  u, v, w, x = U[1,1], U[1,2], U[2,1], U[2,2]
  u2, v2, w2, x2 = u^2, v^2, w^2, x^2
  uv, uw, vx, wx = u * v, u * w, v * x, w * x
  A[1,1] = (u2 - w2 -v2 + x2) // 2
  A[1,2] = uv - wx
  A[1,3] = (u2 - w2 + v2 - x2) // 2
  A[2,1] = uw - vx
  A[2,2] = (u * x + v * w)
  A[2,3] = uw + vx
  A[3,1] = (u2 + w2 - v2 - x2) // 2
  A[3,2] = uv + wx
  A[3,3] = (u2 + w2 + v2 + x2) // 2
  A
end

# transform generating matrix from 3x3 to 2x2 form
function tabcmn_genmat{T<:Integer}(A::AbstractArray{T,2})
  size(A) = (3, 3) || throw(ArgumentError("3x3 matrix required"))
  U = similar(A, (2,2))
  a11, a13, a31, a33 = A[1,1], A[1,3], A[3,1], A[3,3]

  u2 =  a11 + a13 + a31 + a33
  w2 = -a11 - a13 + a31 + a33
  v2 = -a11 + a13 - a31 + a33
  x2 =  a11 - a13 - a31 + a33
  
  isodd(u2) || isodd(v2) || isodd(w2) || isodd(x2) && throw(ArgumentError("parity error"))
  u2, v2, w2, x2 = u2 ÷ 2, v2 ÷ 2, w2 ÷ 2, x2 ÷ 2
  u, v, w, x = isqrt(u2), isqrt(v2), isqrt(w2), isqrt(x2)
  u*u != u2 || v*v != v2 || w*w != w2 || x*x != x2 && throw(ArgumentError("no square roots"))
 
  uv = A[3,2] + A[1,2]
  wx = A[3,2] - A[1,2]
  uw = A[2,3] + A[2,1]
  vx = A[2,3] - A[2,1]
  a  = A[2,2]
  isodd(uv) || isodd(wx) || isodd(uw) || isodd(vx) && throw(ArgumentError("odd sums"))
  uv, wx, uw, vx = uv ÷ 2, wx ÷ 2, uw ÷ 2, vx ÷ 2
 
  uv*uv != u2*v2 || wx*wx != w2*x2 || uw*uw != u2*w2 || vx*vx != v2*x2 && throw(ArgumentError("products"))
  u*x + v*w != abs(a) && abs(u*x-v*w) != abs(a) && throw(ArgumentError("midpoint dot product") )
  sign(uv) * sign(wx) * sign(uw) * sign(vx) >= 0 || throw(ArgumentError("sign parity"))
  ux = a
  vw = a
  if uv == 0 && uw == 0 && vx == 0 && wx == 0
    # if u != 0 => w == 0 && v == 0 && x maybe != 0
    # then u * x = a[
    if u != 0
	  x = copxsign(x, a)
	elseif v != 0
	  w = copysign(w, a)
	elseif w != 0
	  v = copysign(v, a)
	elseif x != 0
	  u = copysign(u, a)
	end
  elseif u == 0
	v = copysign(v, vx)
    w = copysign(w, wx)
  elseif v == 0
	u = copysign(u, uw)
    x = copysign(x, wx)
  elseif w == 0
	u = copysign(u, uv)
    x = copysign(x, vx)
  elseif x == 0
	v = copysign(v, uv)
    w = copysign(w, uw)
  else
    v = copysign(v, uv)
	w = copysign(w, uw)
	x = copysign(x, w*wx)
  end
   [u v; w x]
end

# standard generating triple in m,n-form
genmat1 = [[2 -1; 1 0], [2 1; 1 0], [1 2; 0 1]]
# alternative generating triple in m,n-form
genmat2 = [[1 1; 0 2], [2 0; 1 -1], [2 0; 1 1]]
# default generators
genmat = genmat1

# integer logarithm to the basis of 3			
function ilog3{T<:Integer}(n::T)
  prod = typeof(n)(3); k = 0
  while prod <= n
    prod *= 3
	k += 1
  end
  mintype(k)(k)
end

# generate full 2x2 - matrix in canonical order (starting with 0)
#
function getmn_full{T<:Integer}(n::T)
  dec = decompose(n)
  TT = mintype(BigInt(floor(norm([2;1]) * 2.4142136 ^ length(dec))))
  A = eye(TT, 2, 2)
  for k in decompose(n)
    A = genmat[k+1] * A
  end
  A
end

# generate n-th pythagorean m,n-form in canonical order (starting with 0)
#
function getmn{T<:Integer}(n::T)
  A = getmn_full(n)
  A[1,1], A[2,1]
end

# generate n-th pythagorean triple in abc-form in canonical order (starting with 0)
#
function getabc{T<:Integer}(n::T)
  tmnabc(getmn(n)...)
end

# generate full 3x3 matrix
function getabc_full{T<:Integer}(n::T)
  tmnabc_genmat(getmn_full(n)) 
end

# Determinant of integer matrix (of small size; complexity n!)
function det{T<:AbstractArray{Int64,2}}(A::T)
  n = size(A,1)
  n != size(A,2) && throw(ArgumentError("require square matrix"))
  n <= 0 && return Int128(1)
  sum = Int128(0)
  par = 1
  for k = 1:n
    sum += det([A[2:n,1:k-1] A[2:n,k+1:n]]) * A[1,k] * par
	par = -par
  end
  sum
end

# Adjugate matrix of integer square matrix A. (inv(A) * det(A))
function adj{T<:AbstractArray{Int64,2}}(A::T)
  n = size(A,1)
  n != size(A,2) && throw(ArgumentError("require square matrix"))
  B = Array{Int128,2}(n,n)
  pari = 1
  for i = 1:n
    par = pari
    for k = 1:n
      B[k,i] = det([A[1:i-1,1:k-1] A[1:i-1,k+1:n]; A[i+1:n,1:k-1] A[i+1:n,k+1:n]]) * par
      par = -par
    end
	pari = -pari
  end
  B
end

# Inverse of integer matrix with unit determinant
function inv{T<:AbstractArray{Int64,2}}(A::T)
  n = size(A,1)
  n != size(A,2) && throw(ArgumentError("require square matrix"))
  B = adj(A)
  det = vecdot(A[1,:], B[:,1])
  abs(det) == 1 ? B ÷ det : B // det
end

# decompose integer n
function decompose{T<:Integer}(n::T)
  n <= 0 && return Int8[]
  k = n - 1
  list = Int8[]
  while k > 0
    push!(list, mod(k-1, 3))
	k = (k-1) ÷ 3
  end
  push!(list, 0)
  reverse(list)
end

# find unique decomposition of m,n
# as Product of generating matrices applied to [1; 0]
function decompose(m::Integer, n::Integer)
  list = Int8[]
  A = [m; n]
  while any( A .> [1; 0] )
	k, A = findpos(A)
	k >= 0 && push!(list, k)
  end
  reverse(list)
end

function decompose(a::Integer, b::Integer, c::Integer)
  m, n = tabcmn(a, b, c)
  decompose(m, n)
end

decompose(t::Tuple) = decompose(t...)

function findpos(A)
  for k = 0:2
    B = inv(genmat[k+1]) * A
    if all( B .>= 0 ) && B[1] > B[2]
	  return k, B
    end
  end
  return -1, [0; 0]
end

# inverse of decompose
function ord(n::Integer...)
  dec = decompose(n...)
  ld = length(dec)
  ld == 0 && return 0
  n = BigInt(0)
  for k = 2:ld
    n = n * 3 + dec[k] + 1
  end
  n += dec[1] + 1
  mintype(n)(n)
end

ord(t::Tuple) = ord(t...)

# Height-Excess representation of Pythagorean Triples
# If 2*h == p * q² and p squarefree and h, k ∊ ℤ we have:
# When e = d*k
# and a, b, c = e + h, e + e²/(2*h), e + h + e²/(2h)
# then a, b ,c is a Pythagorean triple (Wade, D. McCulloch)
# and all PT can be represented by k, h ∊ ℤ

using Primes
# for integer h = p * q² determine squarefree p and q
function split(h::Integer)
  p = one(h)
  if h == 0
    return h, p 
  end
  for kv in Primes.factor(h)
    ( kv.second & 1 == 1 ) && (p *= kv.first)
  end
  q = isqrt(h ÷ p)
  p, q
end

# transform k, h - form to a,b,c - form
function tkhabc{T<:Integer, S<:Integer}(k::T, h::S)
  if h == 0
    return 0, k, k
  end
  p, q = split(mintype(h, 2)(h) * 2)
  e = p * q * k
  ee = p * k * k
  e + h, e + ee, e + ee + h
end

# transform a,b,c - form to k, h - form
function tabckh{T<:Integer, S<:Integer, U<:Integer}(a::T, b::S, c::U)
  h = c - b
  if h == 0
    return c, h
  end
  p, q = split(mintype(h, 2)(h) * 2)
  d = p * q
  (a - h) ÷ d, h
end

end # module

py = pythagorean
