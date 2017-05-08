function Split(a::Union{Float64, AbstractArray{Float64,}})
global Split
#SPLIT        Error-free split a=x+y into two parts.
#
#   x, y = Split(a)
#
#On return, x + y == a and both x and y need at most k bits of mantissa.
#In double precision k=26, in single precision k=12 (counting implicit 1 of normalized).
#Input may be a vector or matrix as well.
#
#Follows T.J. Dekker: A floating-point technique for extending the available
#  precision, Numerische Mathematik 18:224-242, 1971.
#Requires 4 flops for scalar input.
#
#
    const BM64  = 0x000000000000001 << 26
    const BM641 = ~(BM64 - 1)
    x = reinterpret(UInt64, a)
	x &= BM641
    x += x & BM64
	x = reinterpret(Float64, x)
	y = a - x
	x, y
end

function Split(a::Union{Float32, AbstractArray{Float32,}})
    const BM32::UInt32  = 1 << 12
    const BM321::UInt32 = BM32 - 1
    x = reinterpret(UInt32, a)
	x -= x & BM321
	x += x & BM32
	x = reinterpret(Float32, x)
	y = a - x
	x, y
end
