module Register

export Registry, RegType, assign, dependencies, needs

import Base: colon, =>, >, !
using Syntax
import Syntax:	SequenceSetOfType, SequenceSetType, ChoiceType, ReferenceType,
				Member, ExtensionMark


typealias Registry Dict{Symbol, ASN1_TYPE}
typealias RegType ASN1_TYPE

function assign(s::Symbol, t::ASN1_TYPE, reg::Registry)
	sn = string(s)
	haskey(reg, s) && failparse("typename $sn must not be redefined")
	sn != "" && isupper(sn[1]) || failparse("typename $sn does not start with uppercase letter")
	reg[s] = t
end

"generate tagged ASN1-TYPE"
colon{T<:ASN1_TYPE}(tag::Tag, x::T) = x(tag)
colon{T<:ASN1_TYPE,S<:Integer}(tag::Vector{S}, x::T) = length(tag) == 1 && x(Tag(tag[1]))
"generate enumerated from symbol and integer"
colon(s::Symbol, x::Integer) = Enumerated(s, x)

"generate Member form Symbol and ASN1_TYPE"
(=>){T<:ASN1_TYPE}(s::Symbol, x::T) = Member(s, x)
"generate new Member form Symbol and existing Member"
(=>){T<:ASN1_TYPE}(s::Symbol, x::Member{T}) = Member(s, x.member, x.optional)
"generate Member from Member and Option"
(>){T<:ASN1_TYPE}(x::T, opt::Option) = Member(Symbol(), x, opt)
"convert a Symbol into a ReferenceType"
(!)(s::Symbol) = ReferenceType(s)

typealias SyS Set{Symbol}
const EMPTY = SyS()

"find all type symbols, a given type depends on"
dependencies(typ::ASN1_TYPE, reg::Registry)::SyS = dependencies(typ, reg, EMPTY)
"find all type symbols, a given type symbol depends on"
dependencies(ref::Symbol, reg::Registry)::SyS = dependencies(reg[ref], reg, EMPTY)

"find all type symbols, a given type depends on, if not contained in known set"
dependencies(typ::ASN1_TYPE, reg::Registry, known::SyS)::SyS = EMPTY
dependencies(typ::SequenceSetOfType, reg::Registry, known::SyS)::SyS = dependencies(typ.comp, reg, known)
dependencies(typ::SequenceSetType, reg::Registry, known::SyS)::SyS =
	funion(dep_member, typ.comp, reg, known)
dependencies(typ::ChoiceType, reg::Registry, known::SyS)::SyS =
	funion(dep_member, typ.comp, reg, known)

function dependencies(typ::ReferenceType, reg::Registry, known::SyS)::SyS
	s = typ.ref
	if !in(s, known)
		haskey(reg, s) || failparse("missing definition for type $s")
		t = reg[s]
		set = Set([s])
		union(dependencies(t, reg, union(known, set)), set)
	else
		EMPTY
	end
end

dep_member(m::Member, reg::Registry, known::SyS)::SyS = dependencies(m.member, reg, known)
dep_member(m::ExtensionMark, reg::Registry, known::SyS) = EMPTY

"find all type symbols, a given type needs in any case"
needs(typ::ASN1_TYPE, reg::Registry)::SyS = needs(typ, reg, EMPTY)
"find all type symbols, a given type symbol needs in any case"
needs(ref::Symbol, reg::Registry)::SyS = needs(reg[ref], reg, EMPTY)

"find all type symbols, a given type depends on, if not contained in known set"
needs(typ::ASN1_TYPE, reg::Registry, known::SyS)::SyS = EMPTY
function needs(typ::SequenceSetOfType, reg::Registry, known::SyS)::SyS
	minsize(typ) > 0 ? needs(typ.comp, reg, known) : EMPTY
end

needs(typ::SequenceSetType, reg::Registry, known::SyS)::SyS = funion(needs_member, typ.comp, reg, known)
needs(typ::ChoiceType, reg::Registry, known::SyS)::SyS = fintersect(needs_member, typ.comp, reg, known)

function needs(typ::ReferenceType, reg::Registry, known::SyS)::SyS
	s = typ.ref
	if !in(s, known)
		haskey(reg, s) || failparse("missing definition for type $s")
		t = reg[s]
		set = Set([s])
		union(needs(t, reg, union(known, set)), set)
	else
		failparse("recursive use of $s")
	end
end

needs_member(m::Member, reg::Registry, known::SyS)::SyS =
	m.optional == OPTIONAL ? EMPTY : needs(m.member, reg, known)
needs_member(m::ExtensionMark, reg::Registry, known::SyS) = EMPTY

function funion(f::Function, list::MemberList, reg::Registry, known::SyS)
	s = [f(m, reg, known) for m in list]
	isempty(s) ? EMPTY : union(s...)
end

function fintersect(f::Function, list::MemberList, reg::Registry, known::SyS)::SyS
	foldl(intersect, [f(m, reg, known) for m in list])
end

end # module Register
