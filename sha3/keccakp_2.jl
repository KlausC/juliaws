
export Keccak_p, Keccak_f, PrintBitString

typealias State Array{UInt64,2}

debug = false
macro ifdebug(ex::Expr)
  :( debug ? $ex : nothing )
end

type PrintBitString{T}
  len::Int
  data::Union{BitArray{1}, AbstractArray{T,1}, Array{T,1}, SubArray{T,1,Array{T,1}}}
end

PrintBitString(S::BitArray{1}; len = length(S)) = PrintBitString{Bool}(length(S), S)
PrintBitString(S::Array{Bool,1}; len = length(S)) = PrintBitString{Bool}(length(S), BitArray(S))
PrintBitString{T<:Base.BitUnsigned64}(S::AbstractArray{T,1}; len::Int = length(S)*sizeof(T)*8 ) =
  PrintBitString(len, S)

function Base.show{T}(io::IO, P::PrintBitString{T})
  b = P.len
  if typeof(P.data) == BitArray{1}
    A = P.data.chunks
  else
    A = P.data
  end
  n = -fld(-b, 8)
  A = reinterpret(UInt8, A)
  for k = 1:n
	  write(io, "$(uppercase(hex(A[k], 2))) ")
  end
end

# Keccak_p
# Generate Keccak permutation function with pararmeters b and nr
# b must be of form 2^l * 25 where l is in 0:6
#
function Keccak_p(b::Integer, nr::Integer)
  w = b ÷ 25
  (b == w * 25 && ( 0 <= w <= 64)) ||
    throw(ArgumentError("b ($b) is not a positive multiple of 25 <= 1600"))

  l = Int64(log2(w)) # throws InexactError if log2(w) is not integera
  maskw = (UInt64(1) << w ) - 1

  function keccakpnr_debug(S)
  	A = convert_to_State(S)
	B = State(5,5)
	  println("Input of permutation:\n$(PrintBitString(S, len=b))")
      println("Same, with lanes as $(w)- bit words")
	  print_state(A)
	for ir = (12 + 2 * l - nr):(12 + 2*l - 1)
	  Rnd!(A, B, ir)
	end
	S1 = typeof(S) == BitArray{1} ? convert_to_BitArray(A) : convert_to_Array(A, eltype(S))
	  println("State after permutation:\n$(PrintBitString(S1, len=b))")
    S1
  end

  function keccakpnr(S)
  	A = convert_to_State(S)
	B = State(5,5)
	for ir = (12 + 2 * l - nr):(12 + 2*l - 1)
	  Rnd!(A, B, ir)
	end
	S1 = typeof(S) == BitArray{1} ? convert_to_BitArray(A) : convert_to_Array(A, eltype(S))
    S1
  end

  # IMPLEMENTATION FUNCTIONS

function convert_to_State(S::Array{UInt64,1})
  b <= length(S) * 64 || throw(ArgumentError("bit array size must be $b"))

  A = Array{UInt64,2}(5,5)
  k, j = 0, 0
  sk::UInt64 = 0
  for y = 1:5, x = 1:5
    if j == 0
      k += 1
      sk = S[k]
    else
      sk >>= w
    end
    j = mod(j + w, 64)
	A[y,x] = sk & maskw
  end
  A
end

function convert_to_UInt64Array(A::State)
  S = zeros(UInt64, (b+63)÷64)
  k, j = 0, 0
  sk::UInt64 = 0
  @inbounds for y = 1:5, x = 1:5
    if j == 0
      sk = A[y,x]
    else
      sk |= (A[y,x] << j)
    end
    j = mod(j + w, 64)
    if j == 0
      k += 1
	  S[k] = sk
    end
  end
  if j != 0
    S[k+1] = sk
  end
  S
end

function convert_to_State(S::BitArray{1})
  b == length(S) || throw(ArgumentError("bit array size must be $b"))
  convert_to_State(S.chunks)
end

function convert_to_BitArray(A::State)
  S = BitArray(b)
  S.chunks = convert_to_UInt64Array(A)
  S
end

function convert_to_State{T<:Base.BitUnsigned}(S::Array{T,1})
	d
  b <= length(S) * sizeof(T) * 8 || throw(ArgumentError("bit array size must be $b"))
  S1 = zeros(UInt64, (b + 63) ÷ 64)
  copy!(reinterpret(eltype(S), S1), S)
  convert_to_State(reinterpret(UInt64, S1))
end

function convert_to_Array(A::State, T::Type = UInt8)
  n = sizeof(T) * 8
  n = (b + n - 1) ÷ n
  S = Array{T,1}(n)
  copy!(S, 1, reinterpret(T, convert_to_UInt64Array(A)), 1, n)
end

function rotl(u::UInt64, offset::Int)
  offset = mod(offset, w)
  ((u << offset) | (u >> (w-offset))) & maskw
end

# Algorithm 1: θ(A)
function θ!(A::State, As::State)
#	copy!(As, A); return
  C = map( k -> ($)(A[:,k]...), 1:5)
  Cs = map( x -> rotl(x, 1), C)
  D = map( k -> C[mod(k-1, 5)+1] $ Cs[mod(k+1, 5)+1], 0:4)
  for k = 1:5
	As[:,k] = A[:,k] $ D[k]
  end
  nothing
end

# Algorithm 2: ρ(A)
function ρ!(A::State, As::State)
#	copy!(As, A); return
  for x = 1:5, y = 1:5
    As[x,y] = rotl(A[x,y], RhoOffset[x,y])
  end
  nothing
end

# Algorithm 3: π(A)
function π!(A::State, As::State)
#	copy!(As, A); return
  for x = 0:4
    for y = 0:4
      As[y+1,x+1] = A[x+1,mod(3*y+x, 5)+1]
    end
  end
  nothing
end

# Algorithm 4: χ(A)
function χ!(A::State, As::State)
#	copy!(As, A); return
  for x = 0:4
    for y = 0:4
      As[y+1,x+1] = A[y+1,x+1] $ (~(A[y+1,mod(x+1,5)+1]) & A[y+1,mod(x+2,5)+1])
    end
  end
  nothing
end

# Algorithm 7: ι(A, ir) iota
function ι!(A::State, ir::Int)
#	return
  rci = RC[mod(ir,255)+1] & maskw
  A[1,1] $= rci
  nothing
end

# Round ir of Keccak
function Rnd!(A::State, B::State, ir::Int)
  @ifdebug println("\n--- Round $ir ---\n")
  θ!(A, B)
  @ifdebug begin println("After theta:"); print_state(B); end
  ρ!(B, A)
  @ifdebug begin println("B rho:"); print_state(A); end
  π!(A, B)
  @ifdebug begin println("After pi:"); print_state(B); end
  χ!(B, A)
  @ifdebug begin println("After chi:"); print_state(A); end
  ι!(A, ir)
  @ifdebug begin println("After iota:"); print_state(A); end
  nothing
end

print_state(A::Array{UInt64,2}) = println(PrintState(w, A))

### END OF IMPLEMENTATION FUNCTIONS

  debug ? keccakpnr_debug : keccakpnr
end

# Keccak_f
Keccak_f(b::Int) = Keccak_p(b, 12 + 2 * Int64(log2(b÷25)))

# Printing
type PrintState
  w::Int
  A::State
end

function Base.show(stream::IO, P::PrintState)
  A = P.A
  w = P.w
  for y = 1:5
    for x = 1:5
		write(stream, "$(uppercase(hex(A[y,x], (w+3) ÷ 4))) ")
    end
    ( y < 5 ) && write(stream, "\n")
  end
  nothing
end

# some pre-evaluated constants
RC = zeros(UInt64, 255)
RhoOffset = zeros(Int, 5, 5)

let
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

function init_RC()
  for ir = 0:254
    rci = zero(UInt64)
    tl = 1
    for j = 0:6
      if rc(j + 7 * ir)
		  rci |= (UInt64(1) << (tl-1))
      end
      tl += tl
    end
    @ifdebug println("RC[$(dec(ir,2))] = 0x$(hex(rci,16))")
    RC[ir+1] = rci
  end
end

function init_RhoOffset()
  RhoOffset[1,1] = 0
  x, y = 1, 0
  offset = 0
  for t = 1:24
    offset += t
    RhoOffset[y+1,x+1] = mod(offset, 64) 
    x, y = y, mod(2 * x + 3 * y, 5)
  end
  for y = 0:4, x = 0:4
    @ifdebug println("RhoOffset[$x][$y] = $(RhoOffset[x+1,y+1])")
  end
end

  init_RC()
  init_RhoOffset()

end

