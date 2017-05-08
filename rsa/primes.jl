module NextPrimes

export nextprime, randomprime
export random_superprime, is_superprime, next_superprime
import Base: start, done, next, eltype, length
import Base: fldmod, cld

using Primes

"find the next prime number in range r1:(r2+delta)"
function nextprime{T<:Integer, S<:Integer}(r1::T, delta::S)
	r2 = r1 + delta - 1
	(delta > 0 && r2 >= 2) || return T(0)
	r1 <= 2 && return T(2)
	r1 = ((r1 >> 1) << 1) + 1
	r2 = (((r2+1) >> 1) << 1) - 1
	for p in r1:2:r2
		if (p % 3 != 0 || p <= 3) && (p % 5 != 0 || p <= 5)
			if Primes.isprime(p)
				return p
			end
		end
	end
	T(0)
end

"check if prime number p is superprime: (p ÷ 2) - 1 is also prime"
is_superprime(p::Integer) = Primes.isprime(p) && Primes.isprime((p - 1) ÷ 2)

"find next superprime p: p is prime and (p-1) ÷ 2 - 1 is prime"
function next_superprime{T<:Integer, S<:Integer}(r1::T, delta::S)
	found = false
	r2 = r1 + delta - 1
	p = r1 - 1
	while ! found && p <= r2
		p = nextprime(p+1, r2)
		p == 0 && break
		found = Primes.isprime((p-1) ÷ 2)
	end
	found ? p : T(0)
end

"find random prime number in range r1:(r1+delta)"
function randomprime{T<:Integer, S<:Integer}(r1::T, delta::S)
	r2 = r1 + delta - 1
	p1 = rand(r1:r2)
	p = nextprime(p1, r2 - p1 + 1)
	p == T(0) ? nextprime(r1, p1-r1) : p
end

"find random prime number in range r1:(r1+delta)"
function random_superprime{T<:Integer, S<:Integer}(r1::T, delta::S)
	r2 = r1 + delta - 1
	p1 = rand(r1:r2)
	p = next_superprime(p1, r2 - p1 + 1)
	p == T(0) ? next_superprime(r1, p1-r1) : p
end

"Iterator over range of primes"
type PrimeRange{T<:Integer}
	r1::T
	r2::T
end

start(pr::PrimeRange) = nextprime(pr.r1, pr.r2-pr.r1)
done(pr::PrimeRange, s) = s == 0
next(pr::PrimeRange, s) = (s, nextprime(s+1, pr.r2-s))
eltype{T}(pr::PrimeRange{T}) = T
length(pr::PrimeRange) = _length_approximation(pr) + 10

_pi(x::Integer) = x / log(x)
function _length_approximation{T}(pr::PrimeRange{T})
	r1, r2 = pr.r1, pr.r2
	# p1 = T(ceil(log(r1) * r1 * 0.92929))
	# p2 = T(ceil(log(r2) * r2 * 1.1056))
	x = 1.0 / log(r1)
	#T(ceil((_pi(r2) - _pi(r1) + (r2 - r1) * x *(x + 2*x^2 + 6*x^3 + 24* x^4 + 120 * x^5))))
	(r2 - r1 ) ÷ 2 
end

"Euler's phi function. (The number of relatively prime numbers of n <= n)"
phi(n::Integer) = phi(factor(n))
function phi{T<:Integer}(factors::Dict{T,Int})
	prod([ x[1]^(x[2]-1)*(x[1]-1) for x in factors])
end

"Function lambda as defined in RFC 3447 - amended for the case of multiple exponents"
lambda(n::Integer) = lambda(factor(n))
function lambda{T<:Integer}(factors::Dict{T,Int})
	lcm([ x[1]^(x[2]-1)*(x[1]-1) for x in factors])
end
function lambda(factors)
	lcm(collect(factors) - 1)
end

"function count cycles of p modulo n for small values of n, mark tested positions"
function cycles!(p::Int, n::Int, bitv::BitVector)
	k = 1
	p = p % n
	bitv[p] = true
	pk = (p * p) % n
	pk == 0 && return 0
	pk != 0 && (bitv[pk] = true)
	while pk != p && k <= n
		k += 1
		pk = (pk * p) % n
		pk == 0 && return 0
		bitv[pk] = true
	end
	k
end

"find maximal cycle count"
function maxcycle(n::Int)
	bitv = BitVector(n-1)
	fill!(bitv, false)
	cmax = 0
	for p = 1:n-1
		if ! bitv[p]
			c = cycles!(p, n, bitv)
			if c > cmax
				cmax = c
			end
		end
	end
	cmax
end

"check correspondence of maxcycle with lambda and phi"
function checkcycles(io::IO, r::Range)
	for n in r
		d = factor(n)
		pc = length(d)
		if ( pc == sum(values(d)) > 1 )
			maxc = maxcycle(n)
			lamb = lambda(keys(d))
			phiv  = phi(d)
			println(io, "$n $maxc $lamb $phiv $(phiv ÷ lamb) $(keys(d))")
		end
	end
end

end #module
