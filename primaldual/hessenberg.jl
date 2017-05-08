
if ! isdefined(:UpperHessenberg)
immutable UpperHessenberg{T} <: AbstractArray{T,2}
  sub::AbstractVector{T}
  tri::UpperTriangular{T}
  UpperHessenberg{T}(s::AbstractVector{T},t::UpperTriangular{T}) = new(s, t)
end
UpperHessenberg{T}(s::AbstractVector{T},t::UpperTriangular{T}) = UpperHessenberg{T}(s, t)
UpperHessenberg{T}(x::AbstractArray{T}) = UpperHessenberg{T}(diag(x,-1), UpperTriangular(x))
end

if ! isdefined(:LowerHessenberg)
immutable LowerHessenberg{T} <: AbstractArray{T,2}
  sub::AbstractVector{T}
  tri::LowerTriangular{T}
  LowerHessenberg{T}(s::AbstractVector{T},t::LowerTriangular{T}) = new(s, t)
end
LowerHessenberg{T}(s::AbstractVector{T},t::LowerTriangular{T}) = LowerHessenberg{T}(s, t)
LowerHessenberg{T}(x::AbstractArray{T}) = LowerHessenberg{T}(diag(x,1), LowerTriangular(x))
end

import Base: size, getindex, diag, det, \, transpose, ctranspose, inv, full

size{T<:Number}(h::UpperHessenberg{T}) = size(h.tri)
size{T<:Number}(h::LowerHessenberg{T}) = size(h.tri)

getindex{T<:Number}(h::UpperHessenberg{T}, i::Int, j::Int) = i > j + 1 ? zero(T) :
                                                        i == j + 1 ? h.sub[j] :
                                                        h.tri[i,j]

getindex{T<:Number}(h::LowerHessenberg{T}, i::Int, j::Int) = j > i + 1 ? zero(T) :
                                                        j == i + 1 ? h.sub[i] :
                                                        h.tri[i,j]

function getindex{T<:Number}(h::UpperHessenberg{T}, i::UnitRange, j::UnitRange)
  if  i == j
    UpperHessenberg(h.sub[i.start:i.stop-1], UpperTriangular(h.tri[i,j]))
  else
    h[i,j]
  end
end

function getindex{T<:Number}(h::LowerHessenberg{T}, i::UnitRange, j::UnitRange)
  if  i == j
    LowerHessenberg(h.sub[i.start:i.stop-1], LowerTriangular(h.tri[i,j]))
  else
    h[i,j]
  end
end

ctranspose{T<:Number}(h::UpperHessenberg{T}) = LowerHessenberg(conj(h.sub), h.tri')
ctranspose{T<:Number}(h::LowerHessenberg{T}) = UpperHessenberg(conj(h.sub), h.tri')
transpose{T<:Number}(h::UpperHessenberg{T}) = LowerHessenberg(h.sub, h.tri.')
transpose{T<:Number}(h::LowerHessenberg{T}) = UpperHessenberg(h.sub, h.tri.')

function diag{T<:Number}(h::UpperHessenberg{T}, i::Int)
  n = size(h,1)
  i < -1  ? zeros(T, n+i) :
  i == -1 ? copy(h.sub) :
  diag(h.tri,i)
end

function diag{T<:Number}(h::LowerHessenberg{T}, i::Int)
  n = size(h,1)
  i > 1  ? zeros(T, n-i) :
  i == 1 ? copy(h.sub) :
  diag(h.tri,i)
end

function det{T<:Number}(h::UpperHessenberg{T})
  n::Int = size(h, 2)
  p::T = one(T)
  c::Vector{T} = zeros(T, n-1)
  s::Vector{T} = zeros(T, n-1)
  b::T = zero(T)
  a::T = zero(T)
  r::T = zero(T)

  @inbounds for j = 1:n
    b = h[1,j]
    for k = 1:j-1 
      a, b = b, h[k+1,j]
	  b = -s[k] * a + c[k] * b
    end
	if j < n
	  a, b = b, h[j+1,j]
	  G, r = givens(a, b, 1, 2)
	  c[j] = G.c
	  s[j] = G.s
	  b = r
	end
	p *= b
  end
  p
end

det{T<:Number}(h::LowerHessenberg{T}) = det(h.')

function (\){T<:Number}(hh::UpperHessenberg{T}, aa::AbstractVector{T})
  n = size(hh, 1)
  h = UpperHessenberg(copy(hh.sub), copy(hh.tri))
  a = copy(aa)
  @inbounds for k = 1:n-1
    G, r = givens(h[k,k], h[k+1,k], 1, 2)
	h.tri[k,k] = r
	h.sub[k] = zero(T)
	h.tri[k:k+1, k+1:end] = G * h.tri[k:k+1, k+1:end]
	a[k:k+1] = G * a[k:k+1]
  end
  h.tri \ a
end


function (\){T<:Number}(hh::LowerHessenberg{T}, aa::AbstractVector{T})
  n::Int = size(hh, 1)
  h = LowerHessenberg(copy(hh.sub), copy(hh.tri))
  a = copy(aa)
  @inbounds for k = n-1:-1:1
    G, r = givens(h[k+1,k+1], h[k,k+1], 2, 1)
	h.tri[k+1,k+1] = r
	h.sub[k] = zero(T)
	h.tri[k:k+1,1:k] = G * h.tri[k:k+1, 1:k]
	a[k:k+1] = G * a[k:k+1]
  end
  display(h)
  h.tri \ a
end

function div!{T<:Number}(hs::AbstractVector{T}, ht::UpperTriangular{T}, n0, n, aa::AbstractVector{T}, x::AbstractVector{T})
  n = size(hh,1)
  if n0 == n
    x[n0] = ht[n0,n0] \ aa[n0]
  elseif n - n0 == 1
    det = ht[n0,n0] * ht[n,n] - ht[n0,n] * hs[n0]
	x[n0] = (ht[n,n] * aa[n0] - ht[n,n0] * aa[n]) / det
	x[n]  = (ht[n0,n0] * aa[n] - hs[n0] * aa[n0]) / det
  else
    k = div((n0+n), 2)
	ra = n0:k
	rb = k+1:n
	C = hh.tri[ra,rb]
	av = aa[ra]
	bv = aa[rb]
	ak = hh.sub[k]
	ek = map(x->x==k?one(T):zero(T),1:k)
	ek1 = vec(eye(T,n-k,1))
    div!(hs, ht, k+1, n, ek1, x)
	x[ra] = C * x[rb] * ak
	div!(hs, ht, n0, k, x, x)
	de = -x[k] + one(T)
	x0 = div!(A, (av - C * (div!(B, bv))))
	x = x0 + x1 * (x0[k] / de)
	bv[1] -= x[k] * ak
	y = div!(B, bv)
	[x; y]
  end
end

function divi{T<:Number}(hh::UpperHessenberg{T}, aa::AbstractVector{T})
  n = size(hh,1)
  n == size(aa,1) || error("dimension mismatch")
  div!(hh.sub, hh.tri, 1, n, aa, Vector{T}(n))
end

full{T<:Number}(h::UpperHessenberg{T}) = full(sparse(h))
full{T<:Number}(h::UpperHessenberg{T}) = full(sparse(h))
inv{T<:Number}(h::UpperHessenberg{T}) = inv(full(h))
inv{T<:Number}(h::LowerHessenberg{T}) = inv(full(h))


