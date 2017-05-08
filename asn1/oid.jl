
export ObjectIdentifier
import Base: convert, show, hash, ==

"represent an Object Identifier"
immutable ObjectIdentifier
	data::Vector{Int}
	function ObjectIdentifier(arc1::Integer, arc2::Integer, arc::Integer...)
		0 <= arc1 <= 2 || failarg("arc[1] is $arc1 but not in 0..2")
		n = length(arc)
		arc = collect(arc)
		arc2 >= 0 && all( arc .>= 0) || failarg("arcs must not be negative")
		arc1 >= 2 || arc2 <= 39 || failarg("arc2 must be <= 39 for arc1 = $arc2")
		dat = [arc1 * 40 + arc2; arc]
		new(dat)
	end
	ObjectIdentifier(base::ObjectIdentifier, arcs::Integer...) = new([base.data; collect(arcs)])
	ObjectIdentifier{T<:Integer}(dat::Vector{T}) = new(dat)
end

hash(obj::ObjectIdentifier, h::UInt) = hash(obj.data, h)
(==)(a::ObjectIdentifier, b::ObjectIdentifier) = a.data == b.data

"some predifined names of ISO registration"
const ARC_NAMES = Dict{Tuple, String}(
	(0,) => "itu-t*", # final asterisk indicates name need not be followed by parenthesized number
	(0,0) => "recommendation*",
	(0,0,1) => "a*",
	(0,0,2) => "b*",
	(0,0,3) => "c*",
	(0,0,4) => "d*",
	(0,0,5) => "e*",
	(0,0,6) => "f*",
	(0,0,7) => "g*",
	(0,0,8) => "h*",
	(0,0,9) => "i*",
	(0,0,10) => "j*",
	(0,0,11) => "k*",
	(0,0,12) => "l*",
	(0,0,13) => "m*",
	(0,0,14) => "n*",
	(0,0,15) => "o*",
	(0,0,16) => "p*",
	(0,0,17) => "q*",
	(0,0,18) => "r*",
	(0,0,19) => "s*",
	(0,0,20) => "t*",
	(0,0,21) => "u*",
	(0,0,22) => "v*",
	(0,0,23) => "w*",
	(0,0,24) => "x*",
	(0,0,25) => "y*",
	(0,0,26) => "z*",
	(0,1) => "question*",
	(0,2) => "administration*",
	(0,3) => "network-operator*",
	(0,4) => "identified-organization*",
	(1,) => "iso*",
	(1,0) => "standard*",
	(1,2) => "member-body*",
	(1,2,250) => "f",
	(1,2,250,1) => "type-org",
	(1,2,250,1,16) => "ft",
	(1,2,250,1,16,9) => "asn1-book",
	(1,2,840) => "us",
	(1,2,840,11349) => "rsadsi",
	(1,2,840,11349,1) => "pkcs",
	(1,2,840,11349,1,1) => "pkcs-1",
	(1,2,840,11349,1,1,0) => "modules",
	(1,2,840,11349,1,1,0,1) => "pkcs-1",
	(1,2,840,11349,1,1,1) => "rsaEncryption",
	(1,2,840,11349,1,1,2) => "md2WithRSAEncryption",
	(1,2,840,11349,1,1,4) => "md5WithRSAEncryption",
	(1,2,840,11349,1,1,5) => "sha1WithRSAEncryption",
	(1,2,840,11349,1,1,7) => "id-RSAES-OAEP",
	(1,2,840,11349,1,1,8) => "id-mgf1",
	(1,2,840,11349,1,1,9) => "id-pSpecified",
	(1,2,840,11349,1,1,10) => "id-RSASSA-PSS",
	(1,2,840,11349,1,1,11) => "sha256WithRSAEncryption",
	(1,2,840,11349,1,1,12) => "sha384WithRSAEncryption",
	(1,2,840,11349,1,1,13) => "sha512WithRSAEncryption",
	(1,2,840,11349,1,1,14) => "sha224WithRSAEncryption",
	(1,2,840,11349,1,1,15) => "sha512-224WithRSAEncryption",
	(1,2,840,11349,1,1,16) => "sha512-256WithRSAEncryption",
	(1,2,840,11349,2) => "digestAlgorithm",
	(1,2,840,11349,2,2) => "id-md2",
	(1,2,840,11349,2,5) => "id-md5",
	(1,3) => "identified-organization*",
	(1,3,14) => "oiw",
	(1,3,14,3) => "secsig",
	(1,3,14,3,2) => "algorithms",
	(1,3,14,3,2,26) => "id-sha1",
	(2,) => "joint-iso-itu-t*",
	(2,1) => "asn1",
	(2,1,0) => "specification",
	(2,1,0,0) => "modules",
	(2,1,0,0,0) => "iso10646",
	(2,1,0,1) => "characterStrings",
	(2,1,1) => "base-encoding",
	(2,1,2) => "ber-derived",
	(2,1,2,0) => "canonical-encoding",
	(2,1,2,1) => "distinguished-encoding",
	(2,1,3) => "packed-encoding",
	(2,1,3,0) => "basic",
	(2,1,3,0,0) => "aligned",
	(2,1,3,0,1) => "unaligned",
	(2,1,3,1) => "canonical",
	(2,1,3,1,0) => "aligned",
	(2,1,3,1,1) => "unaligned",
	(2,6) => "mhs-motif",
	(2,6,7) => "edims",
	(2,6,7,11) => "bp",
	(2,6,7,11,1) => "edifact-ISO646",
	(2,6,7,11,2) => "edifact-T61",
	(2,6,7,11,3) => "edifact-octet",
	(2,6,7,11,4) => "ansiX12-ISO646",
	(2,6,7,11,5) => "ansiX12-T61",
	(2,6,7,11,6) => "ansiX12-ebcdic",
	(2,9) => "ms",
	(2,16) => "country",
	(2,16,840) => "us",
	(2,16,840,3) => "organization",
	(2,16,840,3,101) => "gov",
	(2,16,840,3,101,3) => "csor",
	(2,16,840,3,101,3,4) => "nistalgorithm",
	(2,16,840,3,101,3,4,2) => "hashalgs",
	(2,16,840,3,101,3,4,2,1) => "id-sha256",
	(2,16,840,3,101,3,4,2,2) => "id-sha384",
	(2,16,840,3,101,3,4,2,3) => "id-sha512",
	(2,16,840,3,101,3,4,2,4) => "id-sha224",
	(2,16,840,3,101,3,4,2,5) => "id-sha512-224",
	(2,16,840,3,101,3,4,2,6) => "id-sha512-256",
) 

function tuple2name(t::Integer...)
	n = length(t)
	num = t[end]
	s = get(ARC_NAMES, t, "")
	length(s) == 0 ? "$num" : s[end] != '*' ? s * "($num)" : s[1:end-1]
end

function show(io::IO, obj::ObjectIdentifier)
	dat = obj.data
	arc1 = dat[1] ÷ 40
	arc2 = dat[1] - arc1 * 40
	n = length(dat)
	names = [ tuple2name(arc1, arc2, dat[2:k]...) for k = 1:n ]
	names = foldl((a::String, b::String) -> a * " " * b, names)
	write(io, "$(tuple2name(arc1)) $names")
	nothing
end

