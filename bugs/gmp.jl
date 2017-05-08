# This file is a part of Julia. License is MIT: http://julialang.org/license

module GMPX


import Base: base, hex

bin1(n::BigInt, p::Int) = base( 2, n, p)
oct1(n::BigInt, p::Int) = base( 8, n, p)
dec1(n::BigInt, p::Int) = base(10, n, p)
hex1(n::BigInt, p::Int) = base(16, n, p)

function base1(b::Integer, n::BigInt, pad::Int)
    nd = ndigits(n, b)
	nd >= pad && return base(b, n)
	corr = BigInt(b) ^ (pad - 1)
	if n.size >= 0
		n += corr
	else
		n -= corr
	end
	s = base(b, n)
	s.data[ifelse(n.size >= 0, 1, 2)] = '0'
	s
end

function base2(b::Integer, n::BigInt, pad::Int)
	2 <= b <= 62 || throw(ArgumentError("base must be 2 ≤ base ≤ 62, got $b"))
	fmt = b == 16 ? "%0*Zx" : b == 10 ? "%0*Zd" : base == 8 ? "%0*Zo" : ""
	pad = max( pad, ndigits(n, b))
	l = n >= 0 ? pad : pad + 1
	p = Vector{UInt8}(l+1) # extra space for \0 and minus-sign
	ccall((:__gmp_snprintf,:libgmp),
		   Ptr{UInt8}, (Ptr{UInt8}, Cint, Ptr{UInt8}, Cint, Ptr{BigInt}), p, l+1, fmt.data, l, &n)
	l = findfirst(p, 0x00) - 1
	String(p[1:l])
end

function hex2(n::BigInt, pad::Int)
	b = IOBuffer()
	res = hex(n)
	diff = pad - length(res)
	for _ in 1:diff
		write(b, "0")
	end
	write(b, res)
	String(b)
end

end
