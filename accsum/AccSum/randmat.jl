# implement the randm,at-algorithm
# The theory is prepared in the article:
#
# A class of arbitrarily ill-conditioned floating-point matrices Siegfried M. Rump
#
# A reference impelmentation for Matlab/Octave exists in the INTMAT package from the same author.
#
#

# Create a nxn matrix with approximate condtion number co, which is exactly
# representable in the given floating point type
function randmat{S<:AbstractFloat}(n::Integer, cond::S)
    k = 2
    P, Q = search_pq(k, cond)
    nmin = -fld(-ndigits(P,2), 52)
    nmax = count_ones(P)
    
    mat = zeros(S, n, n)
	







    mat
end




# construct matrix according to Siegfried M. Rump's suggestion
function constructmat{S<:Real}(m::Int, sigma::Int64, k::Int = 2, v = (1,0,false), tau::Int64 = 0, F::Type{S} = BigInt)
    n, pi, qi, P, Q, det = construct0(m, sigma, k, v, F)
    d = diagm(ones(n-1))
    mat = zeros(S, 2n, 2n)
    mat[1,1:n] = pi
	mat[1, n+1:2n] = qi * k
	mat[2, 1:n] = qi
	mat[2, n+1:2n] = pi
	mat[3:n+1, 1:n-1] = mat[n+2:2n, n+1:2n-1] = d
	mat[3:n+1, 2:n] -= d * sigma
	mat[n+2:2n, n+2:2n] -= d * sigma
    randomize!(mat, P, tau)
end

function constructinv{S<:Real}(m::Int, sigma::Int64, k::Int = 2, v = (1,0,false), tau::Int64 = 0, F::Type{S} = BigInt)
    n, pi, qi, P, Q, det = construct0(m, sigma, k, v, F)
    s = pots(n, sigma)
    si(s,i) = resize0!(s[n-i+1:n], n) * det
    a, b = construct1(n, sigma, pi, qi, P, Q, k, F)
    mat = zeros(S, 2n, 2n)
    sp = s * P
    sq = - s * Q
    mat[1:n,1] = sp
    mat[n+1:2n,1] = sq
    mat[1:n,2] = sq * k
    mat[n+1:2n,2] = sp
    A = s * a'
    B = s * b'
    for i = 1:n-1
        B[:,i] += si(s, i)
    end
    mat[1:n,3:n+1] = B
    mat[n+1:2n,3:n+1] = A
    mat[1:n,n+2:2n] = A * k
    mat[n+1:2n,n+2:2n] = B
    randomizeinv!(mat, P, tau)
end
    
function randomize!{S<:Real}(mat::Array{S,2}, P::BigInt, tau::Int64)
    if tau == 0 return mat end
    z = randomlist(size(minv,1), P, tau)
    for i = 3:size(mat,1)
        mat[i,:] -= mat[1,:] * z[i,1] + mat[2,:] * z[i,2]
    end
    mat
end

function randomizeinv!{S<:Real}(minv::Array{S,2}, P::BigInt, tau::Int64)
    if tau == 0 return minv end
    z = randomlist(size(minv,1), P, tau)
    minv[:,1:2] = minv * z
    minv
end

function randomlist(n2::Int, P::BigInt, tau::Int64)
    rng = MersenneTwister(abs(P))
    nrand() = begin r = rand(rng, -tau+1:tau); r <= 0 ? r-1 : r end
    z = eye(S, n2, 2)
    for i = 3:n2
        z[i,1] = nrand()
        z[i,2] = nrand()
    end
    z
end

function construct0{S<:Real}(m::Int, sigma::Int64, k::Int, v, F::Type{S})
    P, Q, det = pell_pair(k, m, (v[1], v[2]))
    k1 = v[3] ? 1 : k
    pi = expand_p(P, k1, sigma, F)
    qi = expand_p(Q, k1, sigma, F)
    np, nq = length(pi), length(qi)
    n = max(np, nq)
    pi = resize0!(pi, n)
    qi = resize0!(qi, n)
    pi = reverse(pi)
	qi = reverse(qi)
    n, pi, qi, P, Q, det
end

function construct1{S<:Real}(n::Int, sigma::Int64, pi, qi, P, Q, k::Int, F::Type{S})
    a = zeros(F, n-1)
    b = copy(a)
    ai = zero(F)
    bi = ai
    for i = 1:n-1
        ai = ai * sigma + Q * pi[i] - P * qi[i]
        bi = bi * sigma + P * pi[i] - Q * qi[i] * k
        a[n-i] = ai
        b[n-i] = -bi
    end
    reverse!(a), reverse!(b)
end

function pots(n::Integer, sigma::Int)
    s = ones(BigInt, n) * sigma
    s[1] = 1
    s = cumprod(s)
    reverse(s)
end

function resize0!{S}(p::Array{S,1}, n::Integer)
    m = length(p)
    p = resize!(p, n)
	if m < n
	    p[m+1:n] = zeros(S, n-m)
    end
    p
end

function pi_to_p{S<:Real}(pi::AbstractVector{S}, sigma::Int64)
    f(x,y) = x * sigma + y
    foldl(f, BigInt(0), pi)
end
	
import Base.start, Base.next, Base.done

type PellIterator
    k::Int
    x::Int64
    y::Int64
end

function start(pi::PellIterator)
    ( BigInt(1), BigInt(0) )
end

done(pi::PellIterator, st) = false 

function next(pi::PellIterator, st)
	x1, y1 = st
    x2, y2 = pi.x * x1 + pi.k * pi.y * y1, pi.x * y1 + pi.y * x1
    st, (x2, y2)
end

# Pell's equation P² - k * Q² = 1 has infinitely many integer solutions, if
# integer k > 0 is not a square number.
# This function generates an array of solutions for k <= 32
# v = 1, 2 controls the start values of P,Q; either 1,0 or 2,1
# in the second case P² - k * Q² = 2 for all generated P,Q.
#
function pell_numbers(k::Int)
    x, y = pell_fundamental(k)
    PellIterator(k, x, y)
end

function pell_pair(k::Int, n::Int, v = (1,0))
    pit = pell_numbers(k)
    P0, Q0 = minpq(v)
    i = 0
    for (P,Q) in pit
        if i >= abs(n)
            Q *= sign(n)
            P, Q = P0 * P + Q0 * Q * k, Q0 * P + P0 * Q 
            return P,Q, det(P, Q, k)
        end
        i += 1
    end
end

det(P, Q, k) = P * P - Q * Q * k

function pell_fundamental(k::Int)
    if k <= size(PellFundamental,2) && k > 0
        PellFundamental[:,k]
    else
        [0,0]
    end
end

const PellFundamental = [0 0; 3 2; 2 1; 0 0; 9 4; 5 2; 8 3; 3 1; 0 0; 19 6;
						 10 3; 7 2; 649 180; 15 4; 4 1; 0 0; 33 8; 17 4; 170 39; 9 2;
						 55 12; 197 42; 24 5; 5 1; 0 0; 51 10; 26 5; 127 24; 9801 1820; 11 2;
						 1520 273; 17 3]'

# For the given number base k find smallest P,Q that exansion size for k,sigma is >= n
function search_pq(k::Int, sigma::Int, n::Int)
    pit = pell_numbers(k)
	for (P,Q) in pit
		if length(expand_p(P, k, sigma, Float64)) >= n &&
		   length(expand_p(Q, k, sigma, Float64)) >= n

		    return P, Q
		end
	end
end

# For the given number base k find smallest P,Q that the condition number is >= cond
function search_pq{R<:AbstractFloat}(k::Int, cond::R)
    pit = pell_numbers(k)
	for (P,Q) in pit
		if condition(P, Q) >= cond
		    return P, Q
		end
	end
end

condition(P::BigInt, Q::BigInt) = Q == 0 ? one(P) : (div((P^2 -1 ), Q) + P) ^ 2

function condition(m::Int, k::Int = 2, v = (1,0))
    P, Q = pell_pair(k, m, v)
	c = setprecision(53) do
	    float(condition(P, Q))
	end
	c
end

#						 
# Decompose number P as sum 0:n pp[i+1] sigma^i
#
function expand_p{S<:Real}(P::BigInt, β::Int, σ::Int, F::Type{S})
    nP = ndigits(P, 2)
	nσ = ndigits(σ, 2)
	n = div((nP - 1), (nσ - 1)) + 1
    pp = zeros(F, n)
	e = one(F)
    i = 0
    while P != 0
        i += 1
        while β > 1
		    p, r = divrem(P, β)
			if r == 0 P = p; e *= β else break end
	    end
        q, r = divrem(P, σ)
		r = Int(r)
        if (β <= 1) || (rem(q, β) != β - 1) || (q < β)
        #if ( r * 2 <= σ ) || (q < β)
            pp[i] = r * e ; P = q
        else
		    pp[i] = (r - σ) * e ; P = q + 1
		end
    end
    resize!(pp, i)
	pp
end

function minpq(P::Integer, Q::Integer)
    PM, QM = one(P), one(Q)
	P, Q = abs(P), abs(Q)
	while Q >= 0 && P >= 0
	    PM, QM = P, Q
		P, Q = P * 3 - Q * 4, Q * 3 - P * 2
	end
	PM, QM
end

minpq{S<:Integer}(v::Tuple{S,S}) = minpq(v[1], v[2])
