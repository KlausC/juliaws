"""
Chinese Remainder Theorem
"""
module CRT

"""
Find the integer x with x = a[i] modulo q[i] for all i and 0 ≤ x < lcm(q)
given n positve divisors q[i] and remainders a[i].
if no solution exists, return -1.
"""
function chinese_remainder{S<:Integer,T<:Integer}(q::Vector{T}, a::Vector{S})
	
	const ERROR = T(-1)
	const n = length(a)
	all( 0 .<= a .< q ) || throw(ArgumentError("need 0 ≤ remainders < divisors"))
	Q = q[1]
	A = mod(a[1], Q)
	for i = 2:n
		"""
		Invariants:
		Q = lcm(q[1:i-1])
		0 ≤ A < Q
		x = A solves equations:
			x = a[k] modulo q[k] for k in 1:i-1 and 0 ≤ x < Q
			x = A modulo Q (which is equivalent to the previous) set of equations.
		In a single loop run, extend the invariants from i-1 to i.
		try to solve the 2 equations for x
		(1)	x = A modulo Q
		(2)	x = a[i] modulo q[i]
		"""
		qi = q[i]
		ai = mod(a[i], qi)	# enforce 0 ≤ x < qi
		# println("round $i:")
		# printstat([A; ai], [Q; qi])
		c = gcd(Q, qi)	# c == 1 in the classical case of all q[i] coprime
		b = 0
		if c != 1
			"""
			transform equations (1) and (2) by subtracting b, then divide by c.
			new varible y replacing x: y = (x - 3) ÷ c --- x = y * c + b
			That is only possible if A = a[i] modulo c
			"""
			A, b = fldmod(A, c)		# find new A and b: old A = new A * c + b
			ai, bb = fldmod(ai, c)
			b == bb || return ERROR	# no solution possible if A ≠ b modulo c
			# b == bb || println("assumed a[$i] ($(a[i])) = $b modulo $c")
			qi ÷= c
			Q ÷= c	# strictly fld, ÷ is ok because c > 0. A < Q still valid.
			# printstat([A; ai], [Q; qi], c, b)
		end
		if A != ai
			h = invmod(Q, qi) * (ai - A)	# obtained by subtracting (1) from (2)
			h = mod(h, qi)					# calculation modulo qi: enforce 0 ≤ h < qi
			A += h * Q		# 0 ≤ new A ≤ old A + (qi-1)*Q < Q + (qi-1)*Q = qi * Q
		end
		if c == 1
			Q *= qi			# new Q = lcm(old Q, qi) because Q nd qi are coprime
		else
			A = c * A + b	# 0 ≤ new A < c * (qi * old Q - 1) + c = new Q
			Q *= qi * c		# new Q = lcm(oldQ, qi) because gcd is divided out
		end
		# println("$A mod $Q = $A")
		# at the end of the loop the invariants hold as proved in the comments
	end
	A
end

"implementation using recursion and splitting data into halves"
function chinese_remainder_r{S<:Integer,T<:Integer}(q::AbstractVector{T}, a::AbstractVector{S})
	all( 0 .<= a .< q ) || throw(ArgumentError("need 0 ≤ remainders < divisors"))
	try
		q = copy(q)
		a = copy(a)
		_chinese_remainder_p(q, a, 1, length(q), 1)
		a[1]
	catch ex
		! isa(ex, T) && rethrow()
		ex
	end
end

#"working horse avoiding allocation of new arrays"
function _chinese_remainder_p{S<:Integer,T<:Integer}(q::AbstractVector{T}, a::AbstractVector{S}, k1::Int, n::Int, s::Int)
	
	const ERROR = T(-1)
	# println("B: $k1 $n $s $q $a")
	if n == 2
		ps = q[k1]
		qs = q[k1+s]
		ds = a[k1+s] - a[k1]
		c = gcd(ps, qs)
		mod(ds, c) == 0 || throw(ERROR)
		if c != 1
			ps ÷= c
			qs ÷= c
			ds ÷= c
		end
		u = invmod(ps, qs)
		a[k1] += q[k1] * mod(u * ds, qs)
		q[k1] *= qs
		nothing
	elseif n > 2
		n1 = prevpow2(n-1) # previous lower power of 2 (n+1)÷2 is another option
		n2 = n - n1
		k2 = n1 * s + k1
		_chinese_remainder_p(q, a, k1, n1, s)
		_chinese_remainder_p(q, a, k2, n2, s)
		_chinese_remainder_p(q, a, k1, 2, k2 - k1)
	end
	# println("E: $k1 $n $s $q $a")
end

function printstat(a, q, c = 1, b = 0)
	var = c == 1 ? "x" : "y"
	println("$var mod $(q[1]) = $(a[1])")
	println("$var mod $(q[2]) = $(a[2])")
	if c != 1
		println("gcd = $c, common remainder = $b")
	end
end

"Chinese remainder calculation according to Knuth 2. "
function chinese_remainder_knuth{T<:Integer, S<:Integer}(coprimes::Vector{S}, remainders::Vector{T})
	# check coprimes are coprime
	n = length(coprimes)
	length(remainders) == n || throw(ArgumentError("arguments need to be same size"))
	M = lcm(coprimes)
	prod(coprimes) == M || throw(ArgumentError("arguments need to be pairwise coprime"))
	
	function summand(i::Int)
		mi = coprimes[i]
		Mm = M ÷ mi
		invmod(Mm, mi) * remainders[i] * Mm
	end
	x = reduce((+), [summand(i) for i = 1:n])
	mod(x, M)
end

"Calculate c := gcd(u, v) and invmod(u÷c, v÷c)"
function gcdinvmod(u::Integer, v::Integer)
	c, x, y = gcdx(u, v)
	c, mod(x, y ÷ c)
end

end # module
