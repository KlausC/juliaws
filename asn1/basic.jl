module ASN1

export asnwriteval

include("tag.jl")
include("oid.jl")
include("context.jl") 

failarg(x) = throw(ArgumentError(x))
failparse(x) = throw(ParseError(x))

const allow_short_length = true

function asnwriteeoc(io::IO)
	write(io, [UInt8(UT_EOC); UInt8(UT_EOC)])
end

function asnreadtail_eoc(io::IO)
	oct = read(io, UInt8)
	oct == 0x00 || failparse("missing second 0 in EOC marker")
end

const DEFAULT_TAGS = [
	BitVector => UT_BIT_STRING,
	AbstractArray{UInt8,1} => UT_OCTET_STRING
   ]

" deliver default tag for standard constructed types"
function default_tag(val::Any)
	for (typ, tag) in DEFAULT_TAGS
		if isa(val, typ)
			return tag
		end
	end
	return UT_EOC
end

"create OctetString from integer in base 128"
function octet128(val::Integer)
	n = (ndigits(val, 2) + 6) ÷ 7
	oa = OctetString(n)
	oa[n] = UInt8(val & 0x7f)
	for k = n-1:-1:1
		val >>= 7
		oa[k] = UInt8(val & 0x7f) | 0x80
	end
	oa
end

"create concatenation of octet128 encoded integers"
function octet128{T<:Integer}(intseq::AbstractVector{T})
	n = length(intseq)
	res = Vector{UInt8}()
	for k = 1:n
		oa = octet128(intseq[k])
		append!(res, oa)
	end
	res
end

"read integer in base 128 format"
function asnread128(io::IO)::Integer
	oct = read(io, UInt8)
	val = oct & 0x7f
	k = 1
	while (oct & 0x80) != 0x00
		k += 1
		oct = read(io, UInt8)
		if sizeof(val) < k
			val = widen(val)
		end
		val = (val << 7) | (oct & 0x7f)
	end
	val
end

"write tag field"
function asnwritetag(io::IO, tag::Tag)
	tag.val >= 0 || failarg("tag value must be >= 0")
	oct = tag.typ
	tagval = tag.val
	if tagval < 31
		c = write(io, UInt8(oct | tagval))
	else
		c = write(io, UInt8(oct | 0x1f))
		oa = octet128(tagval)
		c += write(io, oa)
	end
	c
end

"read tag field"
function asnreadtag(io::IO)::Tag
	oct = read(io, UInt8)
	tt = oct & 0xe0
	val = oct & 0x1f
	if val == 0x1f
		val = asnread128(io)
	end
	Tag(tt & 0xc0, tt & 0x20, val)
end

"write length field"
function asnwritelength(io::IO, val::Integer)
	val >= -1 || failarg("tag value must be >= 0")
	if 0 <= val <= 0x7f
		write(io, UInt8(val))
	elseif val > 0
		n = (ndigits(val, 16) + 1) ÷ 2
		if n < 0x7f && allow_short_length
			oa = OctetString(n+1)
			oa[1] = UInt8(n) | 0x80
			for k = n+1:-1:2
				oa[k] = UInt8(val & 0xff)
				val >>= 8
			end
			write(io, oa)
		else
			write(io, 0x80) # TODO or throw exception? - for future extensions
		end
	else
		write(io, 0x80)
	end
end

"read length field"
function asnreadlength(io::IO)::Integer
	oct = read(io, UInt8)
	val = oct
	(oct < 0x80) && return val
	(oct == 0x80) && return -1
	(oct == 0xff) && return -1 # TODO or throw exception? - for future extensions
	n = oct & 0x7f
	val = 0x00
	for k = 1:n
		oct = read(io, UInt8)
		if sizeof(val) < k
			val = widen(val)
		end
		val = (val << 8) | oct
		k += 1
	end
	val
end

"write Bool"
function asnwriteval(io::IO, val::Bool, tag::Tag = Tag(UT_BOOL))
	c = asnwritetag(io, tag)
	c += asnwritelength(io, 1)
	c += write(io, ifelse(val, 0xff, 0x00))
end

"read tail of Bool - tag already done"
function asnreadtail_bool(io::IO, ::Context)::Bool
	n = asnreadlength(io)
	n == 1 || failparse("reading bool length")
	oct = read(io, UInt8)
	oct != 0x00
end

"write NULL"
function asnwritenull(io::IO, tag::Tag = Tag(UT_NULL))
	c = asnwritetag(io, tag)
	c += asnwritelength(io, 0)
	c
end

"read tail of NULL - tag already done"
function asnreadtail_null(io::IO, ::Context)
	n = asnreadlength(io)
	n == 0 || failparse("reading NULL length")
end

const ONES_COMPLEMENT = (-1) & 1 == 0

"write Integer"
function asnwriteval(io::IO, val::Integer, tag::Tag = Tag(UT_INTEGER))
	sig = sign(val)
	n = (ndigits(val, 16) + 1 ) ÷ 2 
	if ONES_COMPLEMENT
		val += 1 # transform to 2s complement bits
	end
	oa = OctetString(n+1)
	for k = n+1:-1:2
		oa[k] = val & 0xff
		val >>= 8
	end
	oa[1] = 0x00
	s = ifelse((sig >= 0) == (oa[2] >= 0x80), 0, 1)
	c = asnwritetag(io, tag)
	c += asnwritelength(io, n - s + 1)
	c += write(io, view(oa, s+1:n+1))
	c
end

"read tail of Integer - tag already done"
function asnreadtail_int(io::IO, ::Context)::Integer
	n = asnreadlength(io)
	n >= 0 || failparse("reading Integer length")
	n == 0 && return 0
	val = read(io, UInt8)
	sig = val
	mask = ifelse(sig >= 0x80, 0xff, 0x00)
	val $= mask
	k2 = ifelse(val == 0, 1, 2)
	for k = 2:n
		if sizeof(val) < k2
			val = widen(val)
		end
		oct = read(io, UInt8) $ mask
		val = (val << 8) | oct
		k2 += 1
	end
	sig < 0x80 && return val
	val = -signed(val)
	val - one(val)
end

"write start indefinite constructed type"
function asnwritestart(io::IO, tag::Tag)
	c = asnwritetag(io, Tag(tag.typ, constructed, tag.val))
	c += asnwritelength(io, -1)
end

"write end-of-content"
function asnwriteendval(io::IO)
	asnwriteeoc()
end

"read constructed bitstring"
function asnreadtailc_bitstring(io::IO, context::Context)
	n = asnreadlength(io)
	pos1 = position(io)
	pos2 = pos1
	res = BitVector(0)
	stop = n == 0
	while !stop
		tag = asnreadtag(io)
		if tag.typ == universal && tag.val == UT_BITSTRING
			append!(res, asnreadtail_bitstring(io, context))
		elseif n < 0 && iseoc(tag)
			stop = asnreadtail_eoc(io)
		else
			failparse("unexpected tag $tag in constructed bitstring")
		end
		pos2 = position(io)
		stop = stop || ( n >= 0 && pos2 >= pos1 + n )
	end
	res
end

"read constructed octet-string"
function asnreadtailc_octetstring(io::IO, context::Context)
	n = asnreadlength(io)
	pos1 = position(io)
	pos2 = pos1
	res = OctetString(0)
	stop = n == 0
	while !stop
		tag = asnreadtag(io)
		if tag.typ == universal && tag.val == UT_OCTETSTRING
			append!(res, asnreadtail_octetstring(io, context))
		elseif n < 0 && iseoc(tag)
			stop = asnreadtail_eoc(io)
		else
			failparse("unexpected tag $tag in constructed bitstring")
		end
		pos2 = position(io)
		stop = stop || ( n >= 0 && pos2 >= pos1 + n )
	end
	res
end

"read constructed sequence"
function asnreadtailc_sequence(io::IO, context::Context)::Vector
	n = asnreadlength(io)
	pos1 = position(io)
	pos2 = pos1
	res = Vector{Any}(0)
	stop = n == 0
	while !stop
		tag = asnreadtag(io)
		if n < 0 && iseoc(tag)
			stop = asnreadtail_eoc(io)
		else
			tailfunction = selecttail(tag, context)			
			part = tailfunction(io, subcontext(context))
			push!(res, part)
			context = nextcontext(context)
		end
		pos2 = position(io)
		stop = stop || ( n >= 0 && pos2 >= pos1 + n)
	end
	res
end


"open new write IOBuffer for definite length constructed types"
function asnopen_subseq(io::IO, tag::Tag)
	iosub = IOBuffer()
	c = asnwritetag(io, Tag(tag.typ, constructed, tag.val))
	iosub, c
end

"close new write IOBuffer and write to channel"
function asnclose_subseq(io::IO, iosub::Base.AbstractIOBuffer)
	data = takebuff_array!(iosub)
	n = length(data)
	c = asnwritelength(io, n)
	c += write(io, data)
	c
end

"write bit string"
function asnwriteval(io::IO, val::BitVector, tag::Tag = Tag(UT_BIT_STRING))
	c = asnwritetag(io, tag)
	nb = length(val)
	n = (nb + 7) ÷ 8
	ni = n * 8 - nb
	c += asnwritelength(io, n+1)
	c += write(io, UInt8(ni))
	oc = 0x00
	for k = 1:nb
		oc = UInt8((oc << 1) | UInt(val[k]))
		if k % 8 == 0 
			c += write(io, oc)
			println("k = $k, oc = $oc")
			oc = 0x00
		end
	end
	if ni > 0
		c += write(io, oc << ni)
	end
	c
end

"read tail of bit string - tag already done"
function asnreadtail_bitstring(io::IO, ::Context)::BitVector
	n = asnreadlength(io) - 1
	ni = read(io, UInt8)
	nb = n * 8 - ni
	res = BitVector(nb)
	j = 1
	for k = 1:n
		oct = read(io, UInt8)
		while j <= min(nb, k*8)
			res[j] = oct & 0x80 != 0x00 
			oct <<= 1
			j += 1
		end
	end
	res
end

" write octet string"
function asnwriteval(io::IO, val::AbstractVector{UInt8}, tag::Tag = Tag(UT_OCTET_STRING))
	c = asnwritetag(io, tag)
	c += asnwritelength(io, length(val))
	c += write(io, val)
	c
end

"read tail of octet string - tag already done"
function asnreadtail_octetstring(io::IO, ::Context)::OctetString
	n = asnreadlength(io)
	res = OctetString(n)
	c = read!(io, res)
	res
end

"write object identifier"
function asnwriteval(io::IO, val::ObjectIdentifier, tag::Tag = Tag(UT_OBJECT_IDENTIFIER))
	tag.typ & constructed == 0 || failarg("object identifier needs tag in primitive form")
	oa = octet128(val.data)
	c = asnwritetag(io, tag)
	c += asnwritelength(io, length(oa))
	c += write(io, oa)
	c
end

"read tail of object identifier - tag already done"
function asnreadtail_object_identifier(io::IO, ::Context)::ObjectIdentifier
	n = asnreadlength(io)
	dat = Vector{Integer}()
	pos1 = position(io)
	pos2 = pos1 + n
	while pos1 < pos2
		push!(dat, asnread128(io))
		pos1 = position(io)
	end
	ObjectIdentifier(dat)
end

"write UTC time"
function asnwriteval(io::IO, val::DateTime, tag::Tag = Tag(UT_UTCTime))
	tag.typ & constructed == 0 || failarg("UTC time needs tag in primitive form")
	c = asnwritetag(io, tag)
	str = Dates.format(val, "yymmddHHMMSS")
	c += asnwritelength(io, length(str.data))
	c += write(io, str.data)
	c
end

"read tail of UTC time - tag already done"
function asnreadtail_time(io::IO, ::Context)::DateTime
	n = asnreadlength(io)
	oa = Vector{UInt8}(n+2)
	readbytes!(io, view(oa, 3:n+2))
	oa[1:2] = [0x32, 0x30]
	DateTime(String(oa), "yyyymmddHHMMSS"[1:n+2])
end

"write string values for 8-bit character sets"
function asnwrite_characters(io::IO, val::String, tag::Tag, tag2::Tag = tag)
	tag.typ & constructed == 0 || failarg("string needs tag in primitive form")
	findnext([UT_PrintableString, UT_T61String, UT_IA5String], tag.val, 1) > 0 ||
		failarg("tag type for strings is restricted to Printable, T61, IA5")
	c = asnwritetag(io, tag)
	c += asnwritelength(io, length(val.data))
	c += write(io, val.data)
	c
end

"read tail of string - tag already done"
function asnread_characters(io::IO)::String
	n = asnreadlength(io)
	oa = Vector{UInt8}(n)
	readbytes!(io, oa)
	String(oa)
end

"define type simply by universal tag"
const TAILS = Dict{Tag, Function}(
	  Tag(UT_OCTET_STRING) => asnreadtailc_octetstring,
)


end # module
