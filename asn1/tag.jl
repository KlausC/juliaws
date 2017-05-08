
export Tag, cTag, aTag, acTag, isconstructed, tagtype, TAG_NULL
import Base: convert, show, hash, ==

typealias OctetString Array{UInt8,1}

"the four tag types bitmask"
const universal, application, context_specific, private = [0x00, 0x40, 0x80, 0xc0]
"the constructed bitmask"
const simple, constructed = [0x00, 0x20]

"supported values of universal tags"
const UT_EOC, UT_BOOL, UT_INTEGER, UT_BIT_STRING, UT_OCTET_STRING, UT_NULL,
UT_OBJECT_IDENTIFIER, UT_ObjectDescriptor, UT_EXTERNAL, UT_REAL, UT_ENUMERATED,
UT_EMBEDDED_PDV, UT_UTF8String, UT_RELATIVE_OID, UT_SEQUENCE, UT_SET,
UT_NumericString, UT_PrintableString, UT_T61String, UT_VideotexString, UT_IA5String,
UT_UTCTime, UT_GeneralizedTime, UT_GraphicString, UT_VisibleString, UT_GeneralString,
UT_UniversalString, UT_CHARACTER_STRING, UT_BMPString =
[0,1,2,3,4,5, 6,7,8,9,10, 11,12,13,16,17, 18,19,20,21,22, 23,24,25,26,27, 28,29,30]

const UT_TeletexString = UT_T61String
const UT_ISO0646String = UT_VisibleString

"represent tag type and tag value"
immutable Tag
	typ::UInt8
	val::Integer
	function Tag(tt::UInt8, cons::UInt8, val::Integer)
		tt & 0x3f == 0 && 0x00 <= (tt>>6) < 4 || failarg("tag type $tt invalid")
		cons & 0xdf == 0x00 || failarg("tag construction indicator $cons invalid")
		val >= 0 || failarg("tag value $val is negative")
		new(tt | cons, val)
	end
end

uTag(unitag::Int) = Tag(universal, simple, unitag)
ucTag(unitag::Int) = Tag(universal, constructed, unitag)
aTag(unitag::Int) = Tag(application, simple, unitag)
acTag(unitag::Int) = Tag(application, constructed, unitag)
Tag(unitag::Int) = Tag(context_specific, simple, unitag)
cTag(unitag::Int) = Tag(context_specific, constructed, unitag)
tagtype(tag::Tag)::UInt8 = tag.typ & 0xc0
isconstructed(tag::Tag)::Bool = tag.typ & constructed == constructed
iseoc(tag::Tag) = tag.typ == UT_EOC && tag.val == 0
hash(tag::Tag, h::UInt) = hash(tag.val, hash(tag.typ, h))
(==)(a::Tag, b::Tag) = a.typ == b.typ && a.val == b.val
function show(io::IO, tag::Tag)
	ttidx = tag.typ >> 6 + 1
	ttname = ["u", "a", "", "p"][ttidx]
	coidx = (tag.typ >> 5) & 1 + 1
	coname = ["", "c"][coidx]
	if ttidx == 1 && coidx == 1 && tag.val == 0
		print(io, "[]")
	else
		print(io, "[$ttname$coname$(tag.val)]")
	end
	nothing
end

const TAG_NULL = uTag(UT_EOC)

