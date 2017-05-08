
import Base.dump
# access the interior representaion of BigFlooat and BigInt values
function dump(x::BigFloat)
	n = (x.prec + 63) >> 6
    data = reverse(unsafe_wrap(Array, x.d, n))
    sig = x.sign < 0 ? '-' : '+'
	exp = UInt64(x.exp)
	sig, exp, data
end

function dump(x::BigInt)
    n = abs(x.size)
    sig = sign(x.size)
	data = reverse(unsafe_wrap(Array, x.d, n))
    sig, data
end
