"""
RSA public and private keys - variants and conversions.
RSA encoding / decoding algorithms
reconstruct p, q from n, e, and e.
reconstruct p, q from ri and di.
"""
module RSAkeys

export RSAkey, RSApublicKey, RSAprivateKey, rsakeys, bicode, bitLength, rsa_n, rsa_d, rsa_e

using NextPrimes, Primes, CRT
import Base: (==), hash, convert, length

abstract RSAkey

"Public RSA Key"
immutable RSApublicKey{T<:Integer} <: RSAkey
	n::T	# n == p * q (* ...)
	en::T	# public key (renamed to "en" to avoid name clash with Base.e)
end

hash(a::RSApublicKey, h::UInt) = hash(a.en, hash(a.n, h))
rsa_n{T}(key::RSApublicKey{T}) = key.n
length(key::RSApublicKey) = ndigits(key.n, 256)
bitLength(key::RSApublicKey) = ndigits(key.n, 2)

abstract RSAprivateKey <: RSAkey

"Private RSA Key version 0"
immutable RSAprivateKey0{T<:Integer} <: RSAprivateKey
	n::T		# n == p * q (* ...)
	d::T		# private key
	en::T		# public key d * e == 1 mod(phi(p,q,...))	
	mult::Int	# multiplicity of prime factors if known or 0
	function RSAprivateKey0(n, d, en, mult)
		if mult == 0
			try
				pq_from_ned(n, en, d)
				mult = 2
			end
		end
		mult = ifelse(mult == 1, 0, mult)
		new(n, d, en, mult)
	end
end

RSAprivateKey0{T<:Integer}(n::T, d::T, en::T, mult) = RSAprivateKey0{T}(n, d, en, 0)
RSAprivateKey{T<:Integer}(n::T, d::T, en::T) = RSAprivateKey0{T}(n, d, en, 0)
rsa_d{T}(key::RSAprivateKey0{T}) = key.d
rsa_e{T<:RSAprivateKey}(key::T) = key.en
hash(a::RSAprivateKey0, h::UInt) = hash(a.en, hash(a.d, hash(a.n, h)))
multiplicity(a::RSAprivateKey0) = a.mult
defect{T}(a::RSAprivateKey0{T}) = a.mult == 2 ? defect(RSAprivateKey1(a)) : T(1)
rsa_n{T}(key::RSAprivateKey0{T}) = key.n
length(key::RSAprivateKey0) = ndigits(key.n, 256)
bitLength(key::RSAprivateKey0) = ndigits(key.n, 2)

"Private RSA key version 1"
immutable RSAprivateKey1{T<:Integer} <: RSAprivateKey
	q::T
	p::T
	dq::T
	dp::T
	qi::T
	en::T
	RSAprivateKey1(q, p, d, en) = new(q, p, d % (q-1), d % (p-1), invmod(q, p), en)
end

RSAprivateKey{T<:Integer}(q::T, p::T, d::T, en::T) = RSAprivateKey1{T}(q, p, d, en)
rsa_d{T}(key::RSAprivateKey1{T}) = d_from_ri_di([key.q; key.p], [key.dq; key.dp])
hash(a::RSAprivateKey1, h::UInt) = hash(RSAprivateKey0(a), h)
multiplicity(a::RSAprivateKey1) = 2
defect{T}(a::RSAprivateKey1{T}) = gcd(a.p-1, a.q-1)
rsa_n{T}(key::RSAprivateKey1{T}) = key.p * key.q
length(key::RSAprivateKey1) = ndigits(key.p * key.q, 256)
bitLength(key::RSAprivateKey1) = ndigits(key.p * key.q, 2)

"Private RSA key version 1 extended RFC3447"
immutable RSAprivateKey1x{T<:Integer} <: RSAprivateKey
	r::Vector{T} # contains q,   p, r_3, r_4 ...
	d::Vector{T} # contains dq, dp, d_3, d_4, ... 
	t::Vector{T} # contains 1,  qi, t_3, t_4, ... 
	en::T
	function RSAprivateKey1x(r, d, en)
		u = length(r)
		u >= 2 || throw( ArgumentError("at least 2 primes required") )
		da = Vector{T}(u)
		ta = Vector{T}(u)
		ta[1] = 0
		pk = r[1]
		da[1] = d % (pk - 1)
		for k = 2:u
			da[k] = d % (r[k] - 1)
			ta[k] = invmod(pk, r[k])
			pk *= r[k]
		end
		new(r, da, ta, en)
	end
end

" construct RSAprivateKey1x if r[i] and d given"
RSAprivateKey{T<:Integer, S<:Integer}(r::Vector{T}, d::S, en::S) = RSAprivateKey1x{T}(r, d, en)
rsa_d{T}(key::RSAprivateKey1x{T}) = d_from_ri_di(key.r, key.d)
hash(a::RSAprivateKey1x, h::UInt) = hash(RSAprivateKey0(a), h)
multiplicity(a::RSAprivateKey1x) = length(a.r)
defect{T}(a::RSAprivateKey1x{T}) = prod(a.r-1) ÷ lcm(a.r-1)
rsa_n{T}(key::RSAprivateKey1x{T}) = prod(key.r)
length(key::RSAprivateKey1x) = ndigits(prod(key.r), 256)
bitLength(key::RSAprivateKey1x) = ndigits(prod(key.r), 2)

"reconstruct d from values stored in private key"
d_from_ri_di(ri, di) = CRT.chinese_remainder(ri - 1, di)

" try reconstructing p and q from n, e, d"
function pq_from_ned(n, en, d)
	isodd(n) || throw(ArgumentError("n needs to be odd"))
	gcd(d, en) == 1 || throw(ArgumentError("d and e need to be coprime"))
	det = d * en - 1
	nh = (n + 1) ÷ 2
	sqn = isqrt(n)
	if n < 0xffffffffffffffff
		fac = factor(n)
		p = collect(keys(fac))[1]
		q = n ÷ p
	else
		p = q = sqn
		if iseven(p)
			p, q = p-1, q+1
		end
	end
	rmax = nh ÷ max(en, d)
	p0 = 0
	q0 = n
	while p * q !=n && p != p0
		2p < p0 && ( p = p0 ÷ 2; q = n ÷ p)
		q > 2q0 && ( q = q0 * 2; p = n ÷ q) 
		p < 3 && break
		p0, q0  = p, q
		n0 = (p0 - 1) * (q0 - 1) ÷ 2
		# println("det = $det n0 = $n0")
		a, r = continued_fraction(det, n0, rmax)
		# println("a = $a r = $r")
		B = fld(det * r, a)				# 2B approximates (p-1)*(q-1) = n + 1 - p - q
		A = nh - B						# 2A estimates p + q
		y = A^2 - n						# p * q == n
		if y > 0
			x = isqrt(y)
			# println("B = $B A = $(typeof(A))($A) x = $x")
			isodd(A) && isodd(x) && ( x -= 1) # round x to enforce p, q odd
			p, q = A - x, A + x
		else
			p, q = 0, 0
		end
		p <= 2 && ( p = 0 )
		# println("p = $p q = $q")
	end
	p * q == n || throw(ArgumentError("cannot find factorization of n"))
	isprime(p) && Primes.isprime(q) && p > 2 || throw(ArgumentError("cannot find prime factorization of n"))
	p, q
end

"find continued fraction approximating x / y with given maximal denominator"
function continued_fraction(x::Integer, y::Integer, rmax::Integer, round_mode::Int = 0)
	x0, y0 = x, y
	p1, p2, q1, q2 = 1, 0, 0, 1
	pa, qa = p1, q1
	sig = siga = 1
	round_mode = sign(round_mode)
	d() = Float32((widen(y0)*pa - widen(x0)*qa) / qa / y0)

	while y != 0
		sig = -sig
		x1, x2 = fldmod(x, y)
		if sig * round_mode < 0
			x1 = x1 + 1
			x2 = x2 - y
		end
		# println("$x//$y = $x1 rem $x2") 
		p1, p2 = p1 * x1 + p2, p1
		q1, q2 = q1 * x1 + q2, q1
		-rmax <= q1 <= rmax || break
		siga = round_mode == 0 ? sig : round_mode
		pa, qa = p1, q1
		# println("  $pa//$qa $siga $(d())") 
		
		x, y = x2 > 0 ? (y, x2) : (-y, -x2)
	end
	siga = ifelse( y == 0, 0, siga)
	qa > 0 ? (pa, qa, siga) : (-pa, -qa, siga)
end

" RSAprivateKey1 => RSAprivateKey1x"
function convert{T}(::Type{RSAprivateKey1x}, x::RSAprivateKey1{T})
	RSAprivateKey([x.q; x.p], rsa_d(x), x.en)
end

" RSAprivateKey1x => RSAprivateKey1"
function convert{T}(::Type{RSAprivateKey1}, x::RSAprivateKey1x{T})
	length(x.r) == 2 || throw(ArgumentError("cannot convert length > 2"))
	RSAprivateKey(x.r[1], x.r[2], rsa_d(x), x.en)
end

" RSAprivateKey1 => RSAprivateKey0"
function convert{T}(::Type{RSAprivateKey0}, x::RSAprivateKey1{T})
	RSAprivateKey0(x.q * x.p, rsa_d(x), x.en, 2)
end

" RSAprivateKey1x => RSAprivateKey0"
function convert{T}(::Type{RSAprivateKey0}, x::RSAprivateKey1x{T})
	RSAprivateKey0(prod(x.r), rsa_d(x), x.en, multiplicity(x) == 2 ? 2 : 1)
end

" RSAprivateKey0 => RSAprivateKey1x"
function convert{T}(::Type{RSAprivateKey1x}, x::RSAprivateKey0{T})
	multiplicity(x) == 2 || throw(ArgumentError("cannot convert key type 0"))
	p, q = pq_from_ned(x.n, x.en, x.d)
	RSAprivateKey([p; q], x.d, x.en)
end

" RSAprivateKey0 => RSAprivateKey1"
function convert{T}(::Type{RSAprivateKey1}, x::RSAprivateKey0{T})
	multiplicity(x) == 2 || throw(ArgumentError("cannot convert key type 0"))
	p, q = pq_from_ned(x.n, x.en, x.d)
	RSAprivateKey(p, q, x.d, x.en)
end

# equality to compare the different types of private key.
# for the private keys, only n, e, and d are compared, not the respresentation. 

(==)(a::RSApublicKey, b::RSApublicKey) = a.n == b.n && a.en == b.en

(==)(a::RSAprivateKey0, b::RSAprivateKey0) = a.n == b.n && a.d == b.d && a.en == b.en

(==){S<:RSAprivateKey, T<:RSAprivateKey}(a::S, b::T) = RSAprivateKey0(a) == RSAprivateKey0(b)

# key generation using ad-hoc methods

"generate public and private RSA keys with 2 or more prime factors"
function rsakeys(bitsize::Int, gap::Int, u::Int)
	const onbi = one(BigInt)
	# generate smallest p first, then second with gap more bits,
	# the last p is sized as to fill up the remaining bits for the intended
	# key size.
	
	b1 = (4 * bitsize ÷ u - 2 * gap * (u - 1)) ÷ 4 
	b1 = max(b1, 2)
	
	R = onbi
	r = Vector{BigInt}(u)
	for i = 1:u
		p1 = i < u ? p1 = (onbi << b1) + 1 : (onbi<<bitsize) ÷ R
		p1 < 3 && ( p1 = BigInt(3))
		q = randomprime(p1, p1)
		while any( [ gcd(r[k], q) != 1 for k = 1:i-1] )
			q = randomprime(p1, p1)
		end
		r[i] = q
		R *= q
	    b1 = b1 + gap
	end

	n = prod(r)
	phi = lcm(r - 1)
	b2 = min(b1, 16)
	en = BigInt(randomprime((1<<b2)+1, 1<<b2))
	d = invmod(en, phi)

	RSApublicKey(n, en), RSAprivateKey(r, d, en)
end

"generate public and private RSA keys with 2 prime factors"
function rsakeys(bitsize::Int, gap::Int)
	pu, prx = rsakeys(bitsize, gap, 2)
	pu, RSAprivateKey1(prx)
end

# encoding and decoding machines

"encode or decode RSA with public key"
function bicode{T<:Integer,S<:Integer}(key::RSApublicKey{T}, m::S)
	powermod(m, key.en, key.n)
end

"encode or decode RSA with private key version 0"
function bicode{T<:Integer,S<:Integer}(key::RSAprivateKey0{T}, m::S)
	powermod(m, key.d, key.n)
end

"encode or decode RSA with private key version 1"
function bicode{T<:Integer,S<:Integer}(key::RSAprivateKey1{T}, m::S)
	# algo see RFC3477 and http://en.wikipedia.org/wiki/RSA_(cryptosystem)
	p, q = key.p, key.q
	m1 = powermod(m, key.dp, p)
	m2 = powermod(m, key.dq, q)
	h = (m1 - m2) * key.qi
	h = mod(h, p)
	h * q + m2
end

"encode or decode RSA with private key version 1x"
function bicode{T<:Integer,S<:Integer}(key::RSAprivateKey1x{T}, m::S)
	# algo see RFC3477 and http://en.wikipedia.org/wiki/RSA_(cryptosystem)
	# the key. arrays r, d, t are populated: 
	# r = [ q,  p, r3...]
	# d = [dq, dp, d3...]
	# t = [ -, qi, t3...]
	u = length(key.r)
	ma = Vector{T}(u)
	for i = 1:u
		ma[i] = powermod(m, key.d[i], key.r[i])
	end
	R = key.r[1]
	M = ma[1]
	for i = 2:u
		h = (ma[i] - M) * key.t[i]
		h = mod(h, key.r[i])
		M += h * R
		R *= key.r[i]
	end
	M
end

end # module
