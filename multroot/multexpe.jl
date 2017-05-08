
include("multpack.jl")
include("multtest.jl")

module MultExperimental

export enkst, sqsum
export sm, sm1, smn, smn1, sv, svref


import MultiPoly
using MultiPoly
import Polynomials
using Polynomials
import MultRoot
using MultRoot


n, k, s, t = generators(MPoly{Float64}, :n, :k, :s, :t)

A = [n (n+1)*n/2; -(n-(k-1)/2)*k -(n*k-(k-1)*(4k+1)/6)*k]

function enkst(A,n,k,s,t)
    enk(a) = evaluate(a, n, k)
	ank = map(enk, A)
	ank \ [s; t]
end

function sqsum(a::AbstractFloat, b::AbstractFloat, xi, yi, k::Int)

	s1 = sum( (a - b*xi[1:k] -yi[1:k])^2 )
	s2 = sum( (a - b*k -yi[1:k])^2 )
	s = (s1+s2) / 2
end

function sqsum(xi, yi, k::Int)
    n = length(xi)
	if n != length(yi) throw(ArgumentError("length mismatch")) end

    ysum = sum(yi)
	yxsum = sum(xi .* yi)
	ab = enkst(A, n, k, ysum, yxsum)
	a, b = ab[1], ab[2]
	sqsum(a, b, xi, yi, k)
end


sm(p::Poly, m::Int64) = sylves(reverse(p.a), reverse((p').a)/degree(p), m)
sm1(p::Poly, m::Int64) = sylves1(reverse(p.a), reverse((p').a)/degree(p), m)

function scalerow(s::Array{Float64,2})
	m = size(s,1)
	b = Float64[1.0 / norm(s[i,:]) for i = 1:m]
	scale(b, s)
end

smn(p::Poly, m::Int64) = scalerow(sm(p, m))
smn1(p::Poly, m::Int64) = scalerow(sm1(p, m))

sv(p::Poly, m::Int) = zminsv(sm(p,m), 1e-10)

sv(p::Poly) = m::Int -> sv(p,m)[1]


function svref(p::Poly, m::Int64)
    A = sm(p, m)
	s0, x0 = zminsv(A, 1e-10)
	f = reverse(p.a)
	g = reverse((p').a)/(length(p.a)-2)
	h, u, v, res0, res, s1 = MultRoot.gcd_refinement(x0, m, f, g, A)
	s0, s1, res0, res
end

svref(p::Poly) = m::Int -> svref(p,m)

function qrref(p::Poly, m::Int64)
    A = smn1(p, m)
	Q, R = qr(A)
	D = abs(diag(R))
	minimum(D)
end

qrref(p::Poly) = m::Int64 -> qrref(p, m)


end

import MultExperimental
using MultExperimental
