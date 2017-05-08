
module Sha3Hash

include("keccakp.jl")

# Algorithm 8: SPONGE[f, pad, r](N, d)
function sponge_gen(b::Int, f::Function, pad::Function, r::Int, amend::BitArray{})

0 < r <= b || throw(ArgumentError("b ($b) must be > r ($r) > 0"))

# these variables are commonly used be the closures below
# as well as b, f, pad, r, amend
S::BitArray{1} = BitArray(b) 
ptr::Int = 0 

function sponge_init!()
	fill!(S, false)
	ptr = 0
end

function sponge_update!(N::BitArray{1})
	nplus = length(N)
	np = 0
	while np < nplus
		nmax = min(r - ptr, nplus - np)
		println("np= $np nmax = $nmax ptr = $ptr nplus = $nplus")
		S[ptr+1:ptr+nmax] $= N[np+1:np+nmax]
		np += nmax
		ptr += nmax
 		if ptr == r
			S = f(S)
			ptr = 0
		end
	end
end

function sponge_digest!(d::Int)
	if length(amend) > 0
		sponge_update!(amend)
		sponge_update!(pad(r, ptr))
		@assert ptr == 0 "ptr must be 0 after updating pad"
	else
		S = f(S)
	end
	Z = BitArray{1}(d)
	fill!(Z, false)
	nr = min(r, d)
	Z = copy!(Z, 1, S, 1, nr)
	while nr < d
		S = f(S)
		dr = min(r, d - nr)
		copy!(Z, nr + 1, S, 1, dr)
		nr += dr
	end
	sponge_init!()
	Z
end

sponge_init!()			

sponge_update!, sponge_digest!
end

# Algorithm 9: pad10*1(x,m)
function pad101(x::Int, m::Int)
  j = mod(-m-2, x)
  P = BitVector(j+2)
  fill!(P, false)
  P[1] = true
  P[end] = true
  P
end

Keccak(c::Int, a = BitArray([]), b::Int = 1600) = sponge_gen(b, Keccak_f(b), pad101, b - c, a)


_dict = Dict(
			 "SHA3-224" => (448, [0;1], 1600, 224),
			 "SHA3-256" => (512, [0;1], 1600, 256),
			 "SHA3-384" => (768, [0;1], 1600, 384),
			 "SHA3-512" => (1024, [0;1], 1600, 512),
			 "SHAKE128" => (256, [1;1;1;1], 1600, 0),
			 "S24HAKE256" => (512, [1;1;1;1], 1600, 0),
			 "RAWSHAKE128" => (256, [1;1], 1600, 0),
			 "RAWSHAKE256" => (512, [1;1], 1600, 0),
			 "TESTRAW" => (0, [], 1600, 1600) 
	)

type Hasher
	update!::Function
	digest!::Function
	function Hasher(name::String)
		name = uppercase(name)
		haskey(_dict, name) || throw(ArgumentError("Invalid hash type $name"))
		ht = _dict[name]
		c, am, b, d = ht
		up, di = sponge_gen(b, Keccak_f(b), pad101, b - c, BitArray(am))
		new(up, d > 0 ? ()->di(d) : di)
	end
end

end # module

