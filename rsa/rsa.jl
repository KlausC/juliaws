"""
This module implements the methods described in RFC3447 / RFC8017.
The RSA algorithm isself and key generation is outsourced ot RSAkeys
"""
module RSA

export RSAES_AOEP_ENCRYPT

using Nettle

typealias Octet UInt8
typealias OctetString Array{Octet,1}

using RSAkeys, NextPrimes, CRT

typealias MessageRepresentative BigInt
typealias CipherRepresentative BigInt
typealias RSApublicKey RSAkeys.RSApublicKey{BigInt}

immutable RSAoperationError <: Exception
	msg::String
end

fail(s::String) = throw(RSAoperationError(s))

hexdump(s::AbstractArray{UInt8}) = foldl(*, "", [hex(x,2) for x in s])

" 4.1 I2OSP convert integer (BigInt) into array of octets (UInt8) - MSB first"
function i2osp(x::MessageRepresentative, k::Int)
	n = x == 0 ? 0 : ndigits(x, 256)
	n <= k || fail("message representative too big")
	os = zeros(UInt8, k)
	while n > 0 
		os[k] = mod(x, Octet)
		n -= 1
		k -= 1
		x ÷= 256
	end
	os
end

" 4.2 OS2IP - convert array of octets (UInt8) - MSB first - to integer (BigInt)"
function os2ip(os::OctetString)
	k, n = 1, length(os)
    x = BigInt(0)
	while k <= n
		x = x * 256 + os[k]
		k += 1
	end
	x
end

" 5.1.1 RSAEP - encode messagetext with public key"
function RSAEP(key::RSApublicKey, m::MessageRepresentative)
	0 <= m < rsa_n(key) || fail("message representative out of range")
	bicode(key, m)
end

" 5.1.2 RSADP - decode ciphertext with private key"
function RSADP(key::RSAprivateKey, c::CipherRepresentative)
	0 <= c < rsa_n(key) || fail("ciphertext representative out of range")
	bicode(key, c)
end

" 5.2.1 RSASP1 - encode signature with private key"
function RSASP1(key::RSAprivateKey, s::MessageRepresentative)
	0 <= s < rsa_n(key) || fail("message representative out of range")
	bicode(key, s)
end

" 5.2.2 RSAVP1 - decode signature with public key"
function RSAVP1(key::RSApublicKey, s::CipherRepresentative)
	0 <= s < rsa_n(key) || fail("signature representative out of range")
	bicode(key, s)
end

" 7.1.1 Encryption operation RSAES-OAEP-ENCRYPT"
function RSAES_AOEP_ENCRYPT(options, key::RSApublicKey, M::OctetString, L::OctetString = OctetString(0))
	k = length(key)
	mLen = length(M)
	hLen = options.hash_length
	# 1.
	length(L) <= options.hash_input_limit || fail("label too long")
	mLen <= k - 2hLen - 2 || fail("message too long")

	# 2a.
	lHash = options.hash(L)
	# 2b.
	dbLen = k - hLen - 1
	DB = zeros(Octet, dbLen)
	copy!(DB, 1, lHash, 1, hLen)
	copy!(DB, k-mLen-hLen-1, [0x01], 1, 1)
	copy!(DB, k-mLen-hLen, M, 1, mLen)
	# DB = vcat(lHash, PS, 0x01, M)
	# println("DB: $(hexdump(DB))")

	# 2c-e
	seed = options.random(hLen)
	dbMask = options.MGF(seed, dbLen)
	maskedDB = DB $ dbMask
	seedMask = options.MGF(maskedDB, hLen)
	maskedSeed = seed $ seedMask
	# println("maskedSeed: $(hexdump(maskedSeed))")
	# println("maskedDB: $(hexdump(maskedDB))")

	# 2f.
	EM = zeros(Octet, k)
	copy!(EM, 2, maskedSeed, 1, hLen)
	copy!(EM, hLen+2, maskedDB, 1, dbLen)
	# EM = vcat(0x00, maskedSeed, maskedDB)
	# println("EM: $(hexdump(EM))")
	
	# 3a-c.
	m = os2ip(EM)
	c = RSAEP(key, m)
	C = i2osp(c, k)

	# 4.
	C
end

" 7.1.2 Decryption operation RSAES-OAEP-DECRYPT"
function RSAES_OAEP_DECRYPT(options, key::RSAprivateKey, C::OctetString, L::OctetString = OctetString(0))
	hLen = options.hash_length
	k = length(key)
	# 1a-c.
	length(L) <= options.hash_input_limit || fail("decryption error")
	length(C) == k || fail("decryption error")
	k >= 2hLen + 2 || fail("decryption error")

	# 2a-c.
	c = os2ip(C)
	m = RSADP(key, c)
	EM = i2osp(m, k)
	# println("EM: $(hexdump(EM))")

	# 3a-g.
	lHash = options.hash(L)
	Y = EM[1]
	maskedSeed = view(EM, 2:hLen+1)
	maskedDB = getindex(EM, hLen+2:k)
	# println("maskedSeed: $(hexdump(maskedSeed))")
	# println("maskedDB: $(hexdump(maskedDB))")
	#Y, maskedSeed, maskedDB = uncat(EM, 1, hLen, k - hLen - 1)
	seedMask = options.MGF(maskedDB, hLen)
	seed = maskedSeed $ seedMask
	dbMask = options.MGF(seed, k - hLen - 1)
	DB = maskedDB $ dbMask
	# println("DB: $(hexdump(DB))")

	lHash2 = view(DB, 1:hLen)
	y = findnext(x->x==0x01, DB, hLen+1)
	y > 0 && all(view(DB, hLen+1:y-1) .== 0x00) || fail("decryption error")
	M = DB[y+1:end]
	#lHash2, PS, Y2, M = uncat(DB, hLen)
	lHash == lHash2 && Y == 0x00 || fail("decryption error")

	# 4.
	M
end
" 7.2.1 Encryption Operation"
function RSAES_PKCS1_V1_5_ENCRYPT(options, key::RSApublicKey, M::OctetString)
	k = length(key)
	mLen = length(M)

	# 1.
	mLen <= k - 11 || fail("message too long")

	# 2a-b.
	PS = options.random(k - mLen - 3)
	PS[find( x-> x == 0, PS)] = 0x51
	EM = zeros(UInt8, k)
	copy!(EM, 2, [0x02], 1, 1)
	copy!(EM, 3, PS, 1, k - mLen - 3)
	copy!(EM, k - mLen + 1, M, 1, mLen)
	# EM = vcat(0x00, 0x02, PS, 0x00, M)

	# 3a-c.
	m = os2ip(EM)
	c = RSAEP(key, m)
	C = i2osp(c, k)

	# 4.
	C
end

" 7.2.2 Decryption operation"
function RSAES_PKCS1_V1_5_DECRYPT(options, key::RSAprivateKey, C::OctetString)
	k = length(key)

	# 1.
	length(C) == k && k >= 11 || fail("decryption error")

	# 2a-c.
	c = os2ip(C)
	m = RSADP(key, c)
	EM = i2osp(m, k)

	# 3.
	Y, Z = EM[1], EM[2]
	y = findnext(x-> x == 0x00, EM, 3)
	y > 0 || fail("decryption error")
	M = EM[y+1:end]
	#Y, Z, PS, N, M = uncat(EM)
	Y == 0x00 && Z == 0x02 && y >= 11 || fail("decryption error")

	# 4.
	M
end

" 8.1.1 Signature generation operation"
function RSASSA_PSS_SIGN(options, key::RSAprivateKey, M::OctetString)
	k = length(key)
	modBits = bitLength(key)

	# 1.
	EM = EMSA_PSS_ENCODE(options, M, modBits - 1)

	# 2a-c.
	m = os2ip(EM)
	s = RSASP1(key, m)
	S = i2osp(s, k)

	# 3.
	S
end

" 8.1.2 Signature verification operation"
function RSASSA_PSS_VERIFY(options, key::RSApublicKey, M::OctetString, S::OctetString)
	k = length(key)
	modBits = bitLength(key)
	emLen = cld((modBits-1), 8)
	
	# 1.
	length(S) == k || fail("invalid signature")

	# 2a-c.
	s = os2ip(S)
	m = RSAVP1(key, s)
	EM = i2osp(m, emLen)

	# 3.
	Result = EMSA_PSS_VERIFY(options, M, EM, modBits - 1)

	# 4.
	Result == consistent || fail("invalid signature")
end

" 8.2.1 Signature generation operation"
function RSASSA_PKCS1_V1_5_SIGN(options, key::RSAprivateKey, M::OctetString)
	k = length(key)

	# 1.
	EM = EMSA_PKCS1_v1_5_ENCODE(options, M, k)

	# 2a-c.
	m = os2ip(EM)
	s = RSASP1(key, m)
	S = i2osp(s, k)

	# 3.
	S
end

" 8.2.2 Signature verification operation"
function RSASSA_PKCS1_V1_5_VERIFY(options, key::RSApublicKey, M::OctetString, S::OctetString)
	k = length(key)

	# 1.
	length(S) == k || fail("invalid signature")

	# 2a-c.
	s = os2ip(S)
	m = RSAVP1(key, s)
	EM = i2osp(m, k)

	# 3.
	EM2 = EMSA_PKCS1_v1_5_ENCODE(options, M, k)

	EM == EM2 || fail("invalid signature")
	true
end

" 9.1.1 Encoding operation - EMSA-PSS-ENCODE"
function EMSA_PSS_ENCODE(options, M::OctetString, emBits::Integer)

	emLen = (emBits + 7) ÷ 8
	bitmask = 0xff >> (8 * emLen - emBits)

	# 1.
	lengthM = length(M)
	lengthM <= options.hash_input_limit || fail("message too long")
	#2.
	mHash = options.hash(M)
	hLen = options.hash_length
	# 3.
	sLen = options.salt_length
	emBits >= (hLen + sLen + 1) * 8 + 1 || fail("emBits too small")
											  # TODO this check before reading M!
	emLen >= hLen + sLen + 2 || fail("encoding error")
	# 4.
	salt = options.random(sLen)
	# 5.
	MS = zeros(Octet, 8 + hLen + sLen)
	copy!(MS, 9, mHash, 1, hLen)
	copy!(MS, 9 + hLen, salt, 1, sLen)
	# MS = cat(zeros(Octet, 8), mHash, salt)
	# 6.
	H = options.hash(MS)
	# println("sign H: $(hexdump(H))")
	#7.
	PS = zeros(Octet, emLen - sLen - hLen - 2)	# length(PS) may be 0
	# 8.
	dbLen = emLen - hLen - 1
	DB = zeros(Octet, dbLen)
	copy!(DB, emLen - sLen - hLen - 1, [0x01], 1, 1)
	copy!(DB, emLen - sLen - hLen, salt, 1, sLen)
	lz = emLen - hLen - sLen - 2
	# println("sign DB: $lz $(hexdump(DB))")
	#DB = vcat(PS, OctetString([0x01]), salt) # length(DB) = emLen - hLen - 1
	# 9.
	dbMask = options.MGF(H, emLen - hLen - 1)
	# 10.
	maskedDB = dbMask $ DB
	# 11.
	maskedDB[1] &= bitmask
	# 12.
	EM = zeros(Octet, emLen)
	copy!(EM, 1, maskedDB, 1, emLen - hLen - 1)
	copy!(EM, emLen - hLen, H, 1, hLen)
	copy!(EM, emLen, [0xbc], 1, 1)
	#EM = vcat(maskedDB, H, OctetString([0xbc]))
	# println("sign EM: $(hexdump(EM))")
	EM
end

const inconsistent = false
const consistent = true

" 9.1.2 Verification operation - EMSA-PSS-Verify"
function EMSA_PSS_VERIFY(options, M::OctetString, EM::OctetString, emBits::Integer)
	emLen = (emBits + 7) ÷ 8
	bitmask = 0xff >> (8 * emLen - emBits)

	# 1.
	lengthM = length(M)
	lengthM <= options.hash_input_limit || fail("inconsistent signature")
	# 2.
	mHash = options.hash(M)
	hLen = options.hash_length
	
	# 3.
	sLen = options.salt_length
	emLen >= hLen + sLen + 2 || fail("inconsistent signature")
	# 4.
	EM[end] == 0xbc || fail("inconsistent signature")
	# println("verify EM: $(hexdump(EM))")
	# 5.
	maskedDB = view(EM, 1:(emLen - hLen - 1))
	H = view(EM, emLen-hLen:emLen-1)
	# println("verify H: $(hexdump(H))")
	# 6.
	maskedDB[1] & ~bitmask == 0x00 || fail("inconsistent signature")
	# 7.
	dbMask = options.MGF(H, emLen - hLen - 1)
	# 8.
	DB = maskedDB $ dbMask
	# 9.
	DB[1] &= bitmask
	# 10.
	lz = emLen - hLen - sLen - 2
	# println("verify DB: $lz $(hexdump(DB))")

	all(view(DB, 1:lz) .== 0x00) && DB[lz+1] == 0x01 || fail("inconsistent signature")
	# 11.
	salt = view(DB, (lz+2):(emLen - hLen - 1))
	# 12.
	MS = zeros(Octet, 8 + hLen + sLen)
	copy!(MS, 9, mHash, 1, hLen)
	copy!(MS, hLen+9, salt, 1, sLen) 
	#MS = vcat(zeros(Octet, 8), mHash, salt)
	# 13.
	HS = options.hash(MS)
	# 14.
	H == HS || fail("inconsistent signature")
	consistent
end

DER_digest_info(options, H) = H

" 9.2 EMSA-PKCS-v1_5-ENCODE
  This encoding method is deterministic and only has an encoding operation."
function EMSA_PKCS1_v1_5_ENCODE(options, M::OctetString, emLen::Integer)

	# 1.
	H = options.hash(M)
	# 2.
	T = DER_digest_info(options, H)
	tLen = length(T)
	# 3.
	emLen >= tLen + 11 || throw(ArgumentError("intended encoded message length too short"))
	# 4.
	PS = OctetString(emLen - tLen - 3); fill!(PS, 0xff)
	# 5.
	EM = zeros(Octet, emLen)
	copy!(EM, 2, [0x01], 1, 1)
	copy!(EM, 3, PS, 1, emLen - tLen - 3)
	copy!(EM, emLen - tLen + 1, T, 1, tLen)
	#EM = vcat([0x00; 0x01], PS, [0x00], T)
	# 6.
	EM
end

#" B.2.1 MGF1 in options.jl"


include("options.jl")

end #module
