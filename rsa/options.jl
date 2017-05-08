"""
	options for RSA
	1. hash function and related info
	2. MGF - mask generation function
	3. random - random UInt8 array of given size
"""


immutable RSAoption
	hash_length::Int
	hash_input_limit::Int
	salt_length::Int
	der::OctetString
	hash::Function
	MGF::Function
	random::Function
	function RSAoption(hashname::String, mgfhashname::String, sLen::Int, seed::Integer)
		h = MGFoption(hashname)
		mopt = MGFoption(mgfhashname)
		mgf(mseed::AbstractArray{Octet,1}, mLen::Int) = MGF1(mopt, mseed, mLen)

		rng = seed < 0 ? RandomDevice() :
			seed == 0  ? MersenneTwister() : MersenneTwister(seed)
		ran(n::Int) = rand(rng, UInt8, n)
		der = DERstring(hashname)
		new(h.hash_length, h.hash_input_limit, sLen, der, h.hash, mgf, ran)
	end
end

immutable MGFoption
	hash_length::Int
	hash_input_limit::Int
	hash::Function
	function MGFoption(hashname::String)
		h = Hasher(hashname)
		hs = digest!(h)
		hale = length(hs)
		hali = 2^61 - 1 # SHA1-value -- should depend on hash type
		function _hashfu(data::OctetString)
			update!(h, data)
			digest!(h)
		end
		new(hale, hali, _hashfu)
	end
end

const DERdict = Dict(
		"MD2"    => "30 20 30 0c 06 08 2a 86 48 86 f7 0d 02 02 05 00 04 10",
		"MD5"    => "30 20 30 0c 06 08 2a 86 48 86 f7 0d 02 05 05 00 04 10",
		"SHA1"   => "30 21 30 09 06 05 2b 0e 03 02 1a 05 00 04 14",
		"SHA224" => "30 2d 30 0d 06 09 60 86 48 01 65 03 04 02 04 05 00 04 1c",
		"SHA256" => "30 31 30 0d 06 09 60 86 48 01 65 03 04 02 01 05 00 04 20",
		"SHA384" => "30 41 30 0d 06 09 60 86 48 01 65 03 04 02 02 05 00 04 30",
		"SHA512" => "30 51 30 0d 06 09 60 86 48 01 65 03 04 02 03 05 00 04 40",
	"SHA512/224" => "30 2d 30 0d 06 09 60 86 48 01 65 03 04 02 05 05 00 04 1c",
	"SHA512/256" => "30 31 30 0d 06 09 60 86 48 01 65 03 04 02 06 05 00 04 20",
	)

"9.2ote 1."
function DERstring(hashname::String)
	str = get(DERdict, hashname, "")
	io = IOBuffer(str)
	n = (length(str)+1) ÷ 3
	oct = OctetString(n)
	str = filter( x->!isspace(x), str)
	hex2bytes(str)
end

function DER_digest_info(options::RSAoption, H::AbstractArray{Octet,1})
	n = length(options.der)
	hLen = length(H)
	res = OctetString(n + hLen)
	copy!(res, options.der)
	copy!(res, n+1, H, 1, hLen)
	res
end


"B.2.1 MGF1 Mask Generation Function based on a hash function"
function MGF1(options::MGFoption, mgfSeed::AbstractArray{Octet,1}, maskLen::Int)
	hLen = options.hash_length
	
	# 1.
	maskLen <= 2^32 || error("mask too long")
	# 2.
	T = OctetString(maskLen)
	# 3a-b.
	n = cld(maskLen, hLen) - 1
	pos1, pos2 = 1, min(hLen, maskLen)
	for counter = 0:n
		C = i2osp(BigInt(counter), 4)
		hs = options.hash(vcat(mgfSeed, C))
		hLen > (pos2 - pos1 + 1) && (hs = view(hs, 1:pos2-pos1+1))
		T[pos1:pos2] = hs
		pos1 += hLen
		pos2 = min(pos2 + hLen, maskLen)
	end
	# 4.
	T
end

