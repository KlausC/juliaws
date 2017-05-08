# Rotate any unsigned integer type by given amount of shifts.
# Rotation length is the bit-size of input variable.
#
function rotate{T<:Unsigned}(a::T, n::Int)
  sh = mod(n, typeof(a).size * 8)
  b = widen(a) << sh
  bx = reinterpret(T, [b])
  bx[1] | bx[2]
end

function rotate(a::UInt128, n::Int)
  sh = mod(n, 128)
  ax = reinterpret(UInt64, [a])
  if sh > 64
    sh -= 64
    ax[1], ax[2] = ax[2], ax[1]
  end
  b1 = widen(ax[1]) << sh
  b2 = widen(ax[2]) << sh
  b1 | (b2 << 64) | (b2 >> 64)
end

# rotate BigInt variable part of size m by shift of n
function rotate(a::Integer, n::Int, m::Int)
  mm = ( BigInt(1) << m ) - 1
  sh = mod(n, m)
  x1 = (x & mm) << sh
  ( x1 & mm ) | ( x1 >> m )
end
