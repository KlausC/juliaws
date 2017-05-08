
module keccak

type State
  l::Int
  w::Int # == 2^l
  b::Int # == w * 25
  a::Array{UInt64,2}
end

function Base.show(stream::IO, A::State)
  
  for y = 1:5
    for x = 1:5
      write(stream, "$(hex(A.a[y,x], A.w÷4)) ")
    end
    ( y < 5 ) && write(stream, "\n")
  end
end

@inline indexxy(x::Int, y::Int) = (mod(y,5) + 1, mod(x, 5) + 1)
@inline indexz(z::Int, w::Int) = mod(z, w) + 1

function Base.getindex(s::State, x::Int, y::Int, z::Int)
  xx, yy = indexxy(x, y)
  zz = indexz(z, s.w)
  s.a[xx,yy] & mask[zz] == 1
end

function Base.getindex(s::State, x::Int, y::Int)
  xx, yy = indexxy(x, y)
  s.a[xx,yy]
end

function Base.setindex!(s::State, v::Bool, x::Int, y::Int, z::Int)
  xx, yy = indexxy(x, y)
  zz = indexz(z, s.w)
  if v
    s.a[xx,yy] |= mask[zz]
  else
    s.a[xx,yy] &= ~mask[zz]
  end
end

function Base.setindex!(s::State, v::Unsigned, x::Int, y::Int)
  xx, yy = indexxy(x, y)
  s.a[xx,yy] = UInt64(v)
end

function State(S::BitArray{1})
  b = length(S)
  w = b ÷ 25
  l = Int64(log2(w))
  ( 1 << l == w && w * 25 == b) || throw(ArgumentError("Invalid array size $b)"))

  A = State(l, w, b, zeros(UInt64, 5, 5))
  k = 1
  for y = 0:4
    for x = 0:4
	  A[x,y] = S.chunks[k]
	  k += 1
    end
  end
  A
end

function Base.BitArray(A::State)
  S = BitArray(A.b)
  k = 1
  for y = 0:4
    for x = 0:4
	  S.chunks[k] = A[x,y]
      k += 1
    end
  end
  S
end


function rotl(u::UInt64, w::Int, offset::Int)
  offset = mod(offset, w)
  ((u << offset) | (u >> (w-offset))) & (mask[w+1] - 1)
end

# Algorithm 1: θ(A)
function θ!(A::State, As::State)
  C = map( k -> ($)(A.a[:,k]...), 1:5)
  Cs = map( x -> rotl(x, A.w, 1), C)
  D = map( k -> C[mod(k-1, 5)+1] $ Cs[mod(k+1, 5)+1], 0:4)
  for k = 1:5
	As.a[:,k] = A.a[:,k] $ D[k]
  end
end

# Algorithm 2: ρ(A)
function ρ!(A::State, As::State)
  for x = 0:4, y = 0:4
    As[x,y] = rotl(A[x,y], A.w, RhoOffset[x+1,y+1])
  end
end

# Algorithm 3: π(A)
function π!(A::State, As::State)
  for x = 0:4
    for y = 0:4
      As[x,y] = A[mod(3*y+x, 5), x]
    end
  end
end

# Algorithm 4: χ(A)
function χ!(A::State, As::State)
  for x = 0:4
    for y = 0:4
      As[x,y] = A[x,y] $ (~(A[mod(x+1,5),y]) & A[mod(x+2,5),y])
    end
  end
end

# Algorithm 6
function rc(t::Int)
  t = mod(t,255)
  M = 0x8e
  R = 0x80
  for i = 1:t
    R8 = R & 1
    R = R >> 1
    if R8 == 1
      R = R $ M
    end
  end
  Bool(R >> 7)
end

# Algorithm 6: ι(A, ir)
function ι!(A::State, As::State, ir::Int)
  As.a = copy(A.a)
  l = A.l
  rci = RC[mod(ir,255)+1]
  As[0,0] $= rci
end

# Round ir of Keccak
function Rnd!(A::State, ir::Int)
  #D println("\n--- Round $ir ---\n")
  B = State(A.l, A.w, A.b, copy(A.a))
  θ!(B, A)
  #D println("After theta:"); println(A)
  ρ!(A, B)
  #D println("After rho:"); println(B)
  π!(B, A)
  #D println("After pi:"); println(A)
  χ!(A, B)
  #D println("After chi:"); println(B)
  ι!(B, A, ir)
  #D println("After iota:"); println(A)
  A
end

# Keccak_p
function Keccak_p(b::Int, nr::Int)
  w = b ÷ 25
  (b == w * 25) || throw(ArgumentError("b ($b) is not a multiple of 25"))
  l = Int64(log2(w))

  function keccakpnr(S::BitArray{1})
    length(S) == b || throw(ArgumentError("bit string size is not b ($b)"))
  	A = State(S)
    #D println("input of permutation, with lanes as $(A.w) bit words"); println(A)
	for ir = (12 + 2 * l - nr):(12 + 2*l - 1)
	  A = Rnd!(A, ir)
	end
    BitArray(A)
  end
  keccakpnr
end

# Keccak_f
Keccak_f(b::Int) = Keccak_p(b, 12 + 2 * Int64(log2(b÷25)))


# some pre-evaluated constants
mask = map(i -> UInt64(1)<<i, 0:64)
RC = zeros(UInt64, 255)
RhoOffset = zeros(Int, 5, 5)

function init_RC()
  for ir = 0:254
    rci = zero(UInt64)
    tl = 1
    for j = 0:6
      if rc(j + 7 * ir)
        rci |= mask[tl]
      end
      tl += tl
    end
    #D println("RC[$(dec(ir,2))] = 0x$(hex(rci,16))")
    RC[ir+1] = rci
  end
end

init_RC()

function init_RhoOffset()
  RhoOffset[1,1] = 0
  x, y = 1, 0
  offset = 0
  for t = 1:24
    offset += t
    RhoOffset[x+1,y+1] = mod(offset, 64) 
    x, y = y, mod(2 * x + 3 * y, 5)
  end
  for y = 0:4, x = 0:4
    #D println("RhoOffset[$x][$y] = $(RhoOffset[x+1,y+1])")
  end
end

init_RhoOffset()

end # module
