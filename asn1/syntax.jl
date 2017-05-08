module Syntax

export BOOLEAN, NULL, INTEGER, ENUMERATED, BIT_STRING, OCTET_STRING, OBJECT_IDENTIFIER,
		ObjectDescriptor, EXTERNAL, REAL, EMBEDDED_PDV,
		UTF8String,
		RELATIVE_OID, SEQUENCE, SEQUENCE_OF, SET, SET_OF, CHOICE, CHARACTER_STRING,
		NumericString, PrintableString, T61String, VideotexString, IA5String, VisibleString,
		UTCTime, GeneralizedTime, GraphicString, UniversalString, CHARACTER_STRING, BMPString
export EXT, OPTIONAL, Member, MemberList, Constraints, Enumerations, ASN1_TYPE
export Option, Enumerated, StringType, LegacyStringType

import Base: show, colon, (=>), (>), (!)

using CopyModify

include("tag.jl")

# type hierarchy for ASN1 syntax elements
abstract ASN1_TYPE
abstract BASE_TYPE <: ASN1_TYPE
abstract CONSTRUCTED_TYPE <: ASN1_TYPE

"Extension Mark for SELECT, SET, CHOICE, ENUMERATION"
abstract MemberOrExt
immutable ExtensionMark <: MemberOrExt
end

const BL = "  "
const EXT = ExtensionMark()
show(io::IO, em::ExtensionMark, indent::String="", indent2::String=indent) = print(io, "$indent...")

@enum Option OPTIONAL=1 NOT_OPTIONAL=0

immutable Member{T<:ASN1_TYPE} <: MemberOrExt
	name::Symbol
	member::T
	optional::Option
end
Member{T<:ASN1_TYPE}(n::Symbol, m::T) = Member{T}(n, m, false)

function show(io::IO, m::Member, indent::String="", indent2::String=indent)
	print(io, "$indent$(string(m.name)) ")
	show(io, m.member, "", indent2)
	if m.optional == OPTIONAL
		print(io, " OPTIONAL")
	end
end

immutable Enumerated
	name::Symbol
	value::Integer
end

immutable Constraint
end

typealias MemberList{T<:MemberOrExt} Vector{T}
typealias Enumerations Vector{Enumerated}
typealias Constraints Vector{Constraint}
const EMPTY_MEMBERS = Vector{Member}()
const EMPTY_ENUMERATIONS = Enumerations()
const EMPTY_CONSTRAINTS = Constraints()

function show(io::IO, mli::MemberList, indent::String="", indent2::String=indent)
	n = length(mli)
	for k = 1:n
		m = mli[k]
		sep = ifelse(k == n, "", ",")
		show(io, m, indent2)
		println(io, sep)
	end
end

function show(io::IO, eli::Enumerations, indent::String="", indent2::String=indent)
	n = length(eli)
	for k = 1:n
		m = eli[k]
		sep = ifelse(k == n, "", ", ")
		print(io, "$(m.name)($(m.value))$sep")
	end
end

# =====================================
# julia objects representing ASN1 Types
# =====================================

immutable PrimitiveType <: BASE_TYPE
	id::Int
	explicitTag::Tag
	constraints::Constraints
	PrimitiveType(id::Int, exTag::Tag, cs::Constraints) = new(id, exTag, cs)
	PrimitiveType(id::Int) = new(id, TAG_NULL, EMPTY_CONSTRAINTS)
end
immutable LegacyStringType <: BASE_TYPE
	id::Int
	explicitTag::Tag
	constraints::Constraints
	LegacyStringType(id::Int, exTag::Tag, cs::Constraints) = new(id, exTag, cs)
	LegacyStringType(id::Int) = new(id, TAG_NULL, EMPTY_CONSTRAINTS)
end
immutable StringType <: BASE_TYPE
	id::Int
	explicitTag::Tag
	constraints::Constraints
	StringType(id::Int, exTag::Tag, cs::Constraints) = new(id, exTag, cs)
	StringType(id::Int) = new(id, TAG_NULL, EMPTY_CONSTRAINTS)
end
show(io::IO, p::ASN1_TYPE, indent::String="", indent2::String=indent) =
	print(io, "$indent$(stag(p.explicitTag))$(typename(p.id))$(sconstraints(p.constraints))")

stag(tag::Tag) = tag == TAG_NULL ? "" : "$tag "
sconstraints(co::Constraints) = isempty(co) ? "" : " $co"

immutable SequenceSetType <: CONSTRUCTED_TYPE
	id::Int
	explicitTag::Tag
	constraints::Constraints
	comp::MemberList
	SequenceSetType(id::Int, exTag::Tag, cs::Constraints, co::MemberList) =
		new(id, exTag, cs, co)
	SequenceSetType(id::Int) =
		new(id, TAG_NULL, EMPTY_CONSTRAINTS, EMPTY_MEMBERS)
end
function show(io::IO, p::SequenceSetType, indent::String="", indent2::String=indent)
	invoke(show, (IO, ASN1_TYPE, String, String), io, p, indent, indent2)
	println(io, " {")
	show(io, p.comp, indent2 * BL)
	print(io, "$indent2}")
end

immutable SequenceSetOfType <: CONSTRUCTED_TYPE
	id::Int
	explicitTag::Tag
	constraints::Constraints
	comp::ASN1_TYPE
	SequenceSetOfType(id::Int, exTag::Tag, cs::Constraints, co::ASN1_TYPE) =
			new(id, exTag, cs, co)
	SequenceSetOfType(id::Int) = new(id, TAG_NULL, EMPTY_CONSTRAINTS, NULL)
end
minsize(typ::SequenceSetOfType) = 1
maxsize(typ::SequenceSetOfType) = typemax(Int)

function show(io::IO, p::SequenceSetOfType, indent::String="", indent2::String=indent)
	invoke(show, (IO, ASN1_TYPE, String, String), io, p, indent, indent2)
	println(io, " {")
	show(io, p.comp, indent2 * BL)
	print(io, "\n$indent2}")
end


immutable ChoiceType <: CONSTRUCTED_TYPE
	id::Int
	constraints::Constraints
	comp::MemberList
	ChoiceType(id::Int, cs::Constraints, co::MemberList) = new(id, cs, co)
	ChoiceType(id::Int) = new(id, EMPTY_CONSTRAINTS, EMPTY_MEMBERS)
end
function show(io::IO, p::ChoiceType, indent::String="", indent2::String=indent)
	println(io, "$indent$(typename(p.id))$(sconstraints(p.constraints)) {")
	show(io, p.comp, indent2 * BL)
	print(io, "$indent2}")
end

immutable EnumeratedType <: BASE_TYPE
	id::Int
	explicitTag::Tag
	constraints::Constraints
	enumerations::Enumerations
	EnumeratedType(id::Int, exTag::Tag, cs::Constraints, en::Enumerations) =
		new(id, exTag, cs, en)
	EnumeratedType(id::Int) = new(id, TAG_NULL, EMPTY_CONSTRAINTS, EMPTY_ENUMERATIONS)
end
function show(io::IO, p::EnumeratedType, indent::String="", indent2::String=indent)
	print(io, "$indent$(stag(p.explicitTag))$(typename(p.id))$(sconstraints(p.constraints)) {")
	show(io, p.enumerations)
	print(io, "}")
end

immutable ReferenceType <: BASE_TYPE
	ref::Symbol
	explicitTag::Tag
	constraints::Constraints
	ReferenceType(ref::Symbol, exTag::Tag, cs::Constraints) = new(ref, exTag, cs)
	ReferenceType(ref::Symbol) = new(ref, TAG_NULL, EMPTY_CONSTRAINTS)
end
show(io::IO, p::ReferenceType, indent::String="", indent2::String=indent) =
	print(io, "$indent$(stag(p.explicitTag))$(string(p.ref))$(sconstraints(p.constraints))")

# For each ASN1 type there is one object representing default of this type
# An amended modification of each object is created by calling the object with arguments
#

# ===========
# Basic Types
# ===========

const BOOLEAN = PrimitiveType(1)
const INTEGER = PrimitiveType(2)
const BIT_STRING = PrimitiveType(3)
const OCTET_STRING = PrimitiveType(4)
const NULL = PrimitiveType(5)
const OBJECT_IDENTIFIER = PrimitiveType(6)
const ObjectDescriptor = PrimitiveType(7)
const EXTERNAL = PrimitiveType(8)
const REAL = PrimitiveType(9)
const ENUMERATED = EnumeratedType(10)
const EMBEDDED_PDV = PrimitiveType(11)
const RELATIVE_OID = PrimitiveType(13)
const UTCTime = PrimitiveType(23)
const GeneralizedTime = PrimitiveType(24)

# ======================
# Character String Types
# ======================

const NumericString = LegacyStringType(18)
const PrintableString = LegacyStringType(19)
const IA5String = LegacyStringType(20)
const T61String = LegacyStringType(21)
const VideotexString = LegacyStringType(22)
const GraphicString = LegacyStringType(25)
const VisibleString = LegacyStringType(26)
const GeneralString = LegacyStringType(27)
const UniversalString = StringType(28)
const CHARACTER_STRING = LegacyStringType(29)
const BMPString = StringType(30)
const UTF8String = StringType(12)

# =================
# Constructed Types
# =================

const SEQUENCE = SequenceSetType(16)
const SEQUENCE_OF = SequenceSetOfType(14)
const SET = SequenceSetType(17)
const SET_OF = SequenceSetOfType(15)
const CHOICE = ChoiceType(31)

function check_names_tags(co::MemberList)
	names = Vector{String}()
	tags = Vector{Tag}()
	append!(names, filter(s->length(s) > 0, [ p[1] for p in co]))
	append!(tags, filter(t-> t != TAG_NULL, [ tag(p[2]) for p in co]))
	allunique(names) || failparse("identifiers must be unique")
	allunique(tags) || failparse("tags must be unique") 
	co
end

# =====================
"name, tag, prototype, specific functions"
immutable Dispatch
	name::Symbol
	prototype::ASN1_TYPE
	implicitTagValue::Int
end

typename(id::Int) = 1 <= id <= length(DISPATCH) ? (DISPATCH[id].name) : :unknown

"ASN1_TYPE specific dispatch table"
const DISPATCH = Vector{Dispatch}([ 
		Dispatch(:BOOLEAN, BOOLEAN, UT_BOOL),
		Dispatch(:INTEGER, INTEGER, UT_INTEGER), 
		Dispatch(:BIT_STRING, BIT_STRING, UT_BIT_STRING),
		Dispatch(:OCTET_STRING, OCTET_STRING, UT_OCTET_STRING),
		Dispatch(:NULL, NULL, UT_NULL),
		Dispatch(:OBJECT_IDENTIFIER, OBJECT_IDENTIFIER, UT_OBJECT_IDENTIFIER),
		Dispatch(:ObjectDescriptor, ObjectDescriptor, UT_ObjectDescriptor),
		Dispatch(:EXTERNAL, EXTERNAL, UT_EXTERNAL),
		Dispatch(:REAL, REAL, UT_REAL),
		Dispatch(:ENUMERATED, ENUMERATED, UT_ENUMERATED),
		Dispatch(:EMBEDDED_PDV, EMBEDDED_PDV, UT_EMBEDDED_PDV),
		Dispatch(:UTF8String, UTF8String, UT_UTF8String),
		Dispatch(:RELATIVE_OID, RELATIVE_OID, UT_RELATIVE_OID),
		Dispatch(:SEQUENCE_OF, SEQUENCE_OF, UT_SEQUENCE),
		Dispatch(:SET_OF, SET_OF, UT_SET),
		Dispatch(:SEQUENCE, SEQUENCE, UT_SEQUENCE),
		Dispatch(:SET, SET, UT_SET),
		Dispatch(:NumericString, NumericString, UT_NumericString),
		Dispatch(:PrintableString, PrintableString, UT_PrintableString),
		Dispatch(:IA5String, IA5String, UT_IA5String),
		Dispatch(:T61String, T61String, UT_T61String),
		Dispatch(:VideotexString, VideotexString, UT_VideotexString),
		Dispatch(:UTCTime, UTCTime, UT_UTCTime),
		Dispatch(:GeneralizedTime, GeneralizedTime, UT_GeneralizedTime),
		Dispatch(:GraphicString, GraphicString, UT_GraphicString),
		Dispatch(:VisibleString, VisibleString, UT_VisibleString),
		Dispatch(:GeneralString, GeneralString, UT_GeneralString),
		Dispatch(:UniversalString, UniversalString, UT_UniversalString),
		Dispatch(:CHARACTER_STRING, CHARACTER_STRING, UT_CHARACTER_STRING),
		Dispatch(:BMPString, BMPString, UT_BMPString),
		Dispatch(:CHOICE, CHOICE, UT_EOC),
	   ])


let err = filter(k -> DISPATCH[k].prototype.id != k, 1:length(DISPATCH))
	isempty(err) || throw(KeyError("DISPATCH with inconsitent ids: $(err)")) 
end

# =====================

function check_names(comp::MemberList)
	names = Vector{Sting}()
	for co in comp
		append!(names, filter(s->length(s) > 0, [ p[1] for p in co]))
	end
	allunique(names) || failparse("identifiers must be unique")
end

function check_ext(extensible::Bool, comp2::MemberList)
	extensible || length(comp2) == 0 || failparse("extension only if extensible flag set")
end

function leaftypes(typ::DataType)
	function leaftypes(typ::DataType, akku::Vector{DataType})
		if isleaftype(typ)
			push!(akku, typ)
		else
			for t in subtypes(typ)
				leaftypes(t, akku)
			end
		end
	end
	akku = Vector{DataType}()
	leaftypes(typ, akku)
	akku
end

#"construct modified object from original"
for T in leaftypes(ASN1_TYPE)
	@eval function (prod::$T)(cs::Constraint)
		copy_modify(prod, Dict(:constraints => c -> vcat(c,cs)))
	end
	if T != ChoiceType
		@eval function (prod::$T)(tag::Tag)
			copy_modify(prod, Dict(:explicitTag => t -> tag))
		end
	end
end
(prod::SequenceSetType)(co::MemberOrExt...) =
	copy_modify(prod, Dict(:comp => c -> vcat(c,collect(MemberOrExt,co))))

(prod::ChoiceType)(co::MemberOrExt...) =
	copy_modify(prod, Dict(:comp => c -> vcat(c, collect(MemberOrExt,co))))

(prod::SequenceSetOfType)(co::ASN1_TYPE) = copy_modify(prod, Dict(:comp => c -> co))
(prod::EnumeratedType)(en::Enumerated...) =
	copy_modify(prod, Dict(:enumerations => c -> collect(Enumerated, en)))

"throw exeption in case of source errors"
failparse(s::String) = throw(ParseError(s))

end # module
