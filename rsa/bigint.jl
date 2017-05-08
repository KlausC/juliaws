
module BigIntExt

import Base: fldmod, cld

"Implementation of fldmod using gmp - about twice faster than standard"
function fldmod(x::BigInt, y::BigInt)
	z1 = BigInt()
	z2 = BigInt()
	ccall((:__gmpz_fdiv_qr, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &z1, &z2, &x, &y)
	z1, z2
end

"Implementation of cldmod using gmp - this function is not defined in Base"
function cldmod(x::BigInt, y::BigInt)
	z1 = BigInt()
	z2 = BigInt()
	ccall((:__gmpz_cdiv_qr, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &z1, &z2, &x, &y)
	z1, z2
end

"Implementation of cld using gmp - about twice faster than standard"
function cld(x::BigInt, y::BigInt)
	z1 = BigInt()
	ccall((:__gmpz_cdiv_q, :libgmp), Void, (Ptr{BigInt}, Ptr{BigInt}, Ptr{BigInt}), &z1, &x, &y)
	z1
end

end #module
