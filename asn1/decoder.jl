"""
Check compatibility og ASN1 type with concrete julia types.
"""
module Encode

using Syntax
using Register

import Syntax: PrimitiveType, SequenceSetType, SequenceSetOfType, ChoiceType, ReferenceType

type WalkingState
	reg::Registry
	before::Function
	after::Function
	fail::Function
end

function is_compatible(ws::WalkingState, asns::Symbol, obj::Any)
	asn = ws.reg[asns] 
	is_compatible_rec(ws, asn, ojb)
end

function is_compatible_rec(ws::WalkingState, asn::ASN1_TYPE, obj::Any)
	dispatch(asn.id)(ws, asn, obj)
end

function is_compatible_rec(ws::WalkingState, asn::ReferenceType, obj::Any)
	is_compatible_rec(reg, asn.ref, obj)
end

# =============================
# One method for each asn1-type
# =============================

function ic_BOOLEAN(ws::WalkingState, asn::PrimitiveType, obj::Any)
	isa(obj, Bool) ||
	( isa(obj, Integer) && (obj == 0 || obj == 1) ) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_INTEGER(ws::WalkingState, asn::PrimitiveType, obj::Any)
	isa(obj, Integer) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_BIT_STRING(ws::WalkingState, asn::PrimitiveType, obj::Any)
	isa(obj, BitVector) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end


function ic_OCTET_STRING(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Vector{UInt8}) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_NULL(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Void) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_OBJECT_IDENTIFIER(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_ObjectDescriptor(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_EXTERNAL(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_REAL(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, AbstractFloat) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_ENUMERATED(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_EMBEDDED_PDV(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_UTF8String(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, AbstractString) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_RELATIVE_OID(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_SEQUENCE_OF(ws::WalkingState, asn::SequenceSetOfType, ojb::Any)
	isa(obj, AbstractVector) || ws.fail(asn.id)
	for i in 1::length(obj)
		is_compatible_rec(ws, asn.comp, obj[i])
	end
	ws.after(ws, asn, obj)
end

function ic_SET_OF(ws::WalkingState, asn::SequenceSetOfType, ojb::Any)
	isa(obj, Base.AbstractSet) || ws.fail(asn.id)
	for oel in obj
		is_compatible_rec(ws, asn.comp, oel)
	end
	ws.after(ws, asn, obj)
end

function getfield_optional(obj::Any, mem::Member)
	try
		getfield(obj, mem.name)
	catch
		mem.optional == OPTIONAL || rethrow()
		nothing
	end
end 

function ic_SEQUENCE(ws::WalkingState, asn::SequenceSetType, ojb::Any)
	isa(obj, Any) || ws.fail(asn.id)
	for i = 1:length(asn.comp)
		mem = asn.comp[i]
		fld = getfield_optional(obj, sym, mem.optional)
		if fld != nothing
			is_compatible_rec(ws, mem.member, oel)
		end
	end
	ws.after(ws, asn, obj)

end

function ic_SET(ws::WalkingState, asn::SequenceSetType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_NumericString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_PrintableString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_IA5String(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_T61String(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_VideotexString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_UTCTime(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_GeneralizedTime(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_GraphicString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_VisibleString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_GeneralString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_UniversalString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, AbstractString) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_CHARACTER_STRING(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, Bool) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_BMPString(ws::WalkingState, asn::PrimitiveType, ojb::Any)
	isa(obj, AbstractString) || ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

function ic_CHOICE(ws::WalkingState, asn::ChoiceType, ojb::Any)
	failed = true
	for asnel in asn.comp
		try
			is_compatible_rec(ws, asnel.member, obj)
			failed = false
			break
		end
	end
	failed && ws.fail(asn.id)
	ws.after(ws, asn, obj)
end

"dispatch table for walking operation"
DISP_COMPAT = Vector{Function}([
		ic_BOOLEAN,
		ic_INTEGER,
		ic_BIT_STRING,
		ic_OCTET_STRING,
		ic_NULL,
		ic_OBJECT_IDENTIFIER,
		ic_ObjectDescriptor,
		ic_EXTERNAL,
		ic_REAL,
		ic_ENUMERATED,
		ic_EMBEDDED_PDV,
		ic_UTF8String,
		ic_RELATIVE_OID,
		ic_SEQUENCE_OF,
		ic_SET_OF,
		ic_SEQUENCE,
		ic_SET,
		ic_NumericString,
		ic_PrintableString,
		ic_IA5String,
		ic_T61String,
		ic_VideotexString,
		ic_UTCTime,
		ic_GeneralizedTime,
		ic_GraphicString,
		ic_VisibleString,
		ic_GeneralString,
		ic_UniversalString,
		ic_CHARACTER_STRING,
		ic_BMPString,
		ic_CHOICE,
		])

function dispatch(id::Int)
	0 < id <= length(DISP_COMPAT) || error("id $id not in DISP_COMPAT range")
	DISP_COMPAT[id]
end

end # module Compat
