" reverse the order of bits in integer bits types"
module Mirror

export mirror

const MIRRORED_BYTE =
UInt8[0x0, 0x8, 0x4, 0xc, 0x2, 0xa, 0x6, 0xe, 0x1, 0x9, 0x5, 0xd, 0x3, 0xb, 0x7, 0xf] 

"mirror bits of a byte"
function mirror(x::UInt8)::UInt8
	h1 = MIRRORED_BYTE[(x & 0x0f) + 1]
	h2 = MIRRORED_BYTE[(x >> 4) + 1]
	h1 << 4 | h2
end

"mirror bits in unsigned integer"
function mirror{T<:Unsigned}(x::T)::T
	n = sizeof(x)
	y = [ntoh(htol(x))]
	z = reinterpret(UInt8, y)
	for k = 1:n
		z[k] = mirror(z[k])
	end
	y[1]
end

"mirror bits of signed integer"
mirror{T<:Signed}(x::T)::T = signed(mirror(unsigned(x)))

end # module
	
