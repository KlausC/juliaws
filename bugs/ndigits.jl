const ndigits_max_mul = Core.sizeof(Int) == 4 ? 69000000 : 290000000000000000

"ndigits0z failing for UInt128(2)^64"
function ndigits0z_orig(n::Unsigned, b::Int)
	d = 0
	if b < 0
		d = ndigits0znb(signed(n), b)
	else
		b == 2  && return (sizeof(n)<<3-leading_zeros(n))
		b == 8  && return div((sizeof(n)<<3)-leading_zeros(n)+2,3)
		b == 16 && return (sizeof(n)<<1)-(leading_zeros(n)>>2)
		b == 10 && return ndigits0z(n)
		while ndigits_max_mul < n
			n = div(n,b)
			d += 1
		end
		m = 1
		while m <= n
			m *= b
			d += 1
		end
	end
	return d
end

"ndigits0z failing fixed"
function ndigits0z_mod1(n::Unsigned, b::Int)
	if b < 0
		d = ndigits0znb(signed(n), b)
	else
		b == 2  && return (sizeof(n)<<3-leading_zeros(n))
		b == 8  && return div((sizeof(n)<<3)-leading_zeros(n)+2,3)
		b == 16 && return (sizeof(n)<<1)-(leading_zeros(n)>>2)
		b == 10 && return ndigits0z(n)
		b == 16 && return sizeof(n)-(leading_zeros(n)>>3)
		l2n = sizeof(n)<<3 - leading_zeros(n)
		l2b = sizeof(b)<<3 - leading_zeros(b)
		d = div(l2n+l2b, l2b+1)
		n = Int(div(n, typeof(n)(b)^d))
		m = 1 
		while m <= n
			m *= b
			d += 1
		end
	end
	return d
end

"ndigits0z failing fixed"
function ndigits0z_mod2(n::Unsigned, b::Int)
	d = 0
	if b < 0
		d = ndigits0znb(signed(n), b)
	else
		b == 10 && return ndigits0z(n)
		if ispow2(b)
			pow = trailing_zeros(b)
			return div((sizeof(n)<<3-leading_zeros(n)+pow-1), pow)
		end

		while ndigits_max_mul < n
			n = div(n,b)
			d += 1
		end
		m = one(n)
		while m <= n
			m *= b
			d += 1
		end
	end
	return d
end
