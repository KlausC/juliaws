
module keccak

include("keccakp.jl")

# Algorithm 8: SPONGE[f, pad, r](N, d)
function sponge_gen(b::Int, f::Function, pad::Function, r::Int)

0 < r < b || throw(ArgumentError("b ($b) must be > r ($r) > 0"))

function spongefpadd(N::BitArray{1}, d::Int)
  
  P = append!(copy(N), pad(r, length(N)))
  n = length(P) ÷ r
  c = b - r
  j = 0
  S = BitArray(zeros(Bool, b))
  for k = 1:n
    S[1:r] $= P[j+1:j+r]
    S = f(S)
  end
  Z = S[1:r]
  while length(Z) < d
    S = f(S)
    Z = append!(Z, S[1:r])
  end
    length(Z) == d ? Z : Z[1:d]
end

spongefpadd
end

# Algorithm 9: pad10*1(x,m)
function pad101(x::Int, m::Int)
  j = mod(-m-2, x)
  P = BitArray(zeros(Bool,j+2))
  P[1] = true
  P[end] = true
  P
end

Keccak(c::Int, b::Int = 1600) = sponge_gen(b, Keccak_f(b), pad101, b - c)

Keccak256 = Keccak(256)
Keccak448 = Keccak(448)
Keccak512 = Keccak(512)
Keccak768 = Keccak(768)
Keccak1024 = Keccak(1024)

SHA3_224(M::BitArray{1}) = Keccak448(append!(M, [false, true]), 224)
SHA3_256(M::BitArray{1}) = Keccak512(append!(M, [false, true]), 256)
SHA3_384(M::BitArray{1}) = Keccak768(append!(M, [false, true]), 384)
SHA3_512(M::BitArray{1}) = Keccak1024(append!(M, [false, true]), 512)

SHAKE128(M::BitArray{1}, d::Int) = Keccak256(append!(M, BitArray(ones(Bool,4))), d)
SHAKE256(M::BitArray{1}, d::Int) = Keccak512(append!(M, BitArray(ones(Bool,4))), d)
RawSHAKE128(M::BitArray{1}, d::Int) = Keccak256(append!(M, BitArray(ones(Bool,2))), d)
RawSHAKE256(M::BitArray{1}, d::Int) = Keccak512(append!(M, BitArray(ones(Bool,2))), d)

end # module

