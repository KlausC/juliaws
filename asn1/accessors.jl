
module Generic

import Base:	copy

# "make a new instance if mutable data type with fields or return identical"
@generated function Base.copy(x::Any)
	nf = nfields(x)
	(isimmutable(x) || isbits(x) || nf == 0) && return :(x) 
	:(copy_mutable(x), $x, $nf)
	ex
end

function Base.copy(x::Any, dic::Dict{Symbol})
	ty = typeof(x)
	nf = nfields(ty)
	isempty(dic) && return copy(x)
	fields = fieldnames(x)
	y = ccall(:jl_new_struct_uninit, Any, (Any,), ty)::ty
	unmod = true
	for i in 1:nf
		orig = isdefined(x, i) ? getfield(x, i) : nothing
		newv = getval(dic, fields[i], orig)::ty.types[i]
		if newv != nothing
			ccall(:jl_set_nth_field, Void, (Any,Csize_t, Any), y, i-1, newv)
			# setfield!(y, i, newv) fails for immutable
		end
		unmod = unmod && isequal(orig, newv)
	end
	unmod && isimmutable(x) ? x : y
end

function copy_mutable(x::Any, ty::DataType, nf::Int)
	y = ccall(:jl_new_struct_uninit, Any, (Any,), ty)
	for i in 1:nf
		if isdefined(x,i)
			setfield!(y, i, getfield(x, i))
		end
	end
end

@inline function getval(dic::Dict, f::Symbol, val)
	if haskey(dic, f)
		v = dic[f]
		isa(v, Function) ? v(val) : v
	else
		val
	end
end

"define getter and setter functions for all fields of a datatype"
function make_accessors(typ::DataType, gprefix::String="", sprefix::String="")
	if typ.abstract
		for t2 in subtypes(typ)
			make_accessors(t2)
		end
	elseif nfields(typ) > 0
		mod = current_module()
		for name in fieldnames(typ)
			gf = Symbol(gprefix, name)
			expr = :($gf(x::$typ)  = getfield(x, $(QuoteNode(name)) ))
			eval(mod, expr)
			if typ.mutable
				sf = Symbol(sprefix, name)
				eval(mod, :($sf(x::$typ, v) = setfield!(x, $(QuoteNode(name)), v) ))
			end
		end
	end
end

type Undef
end
const UNDEF = Undef() 

"define == and hash functions for a data type - using all fields"
function make_equal{typ}(::Type{typ})

	isleaftype(typ) || return nothing

	mod = current_module()
	fields = fieldnames(typ)

	get = (x::typ, name::Symbol) ->
	try
		eval(mod, name)(x)
	catch
		try
			getfield(x, name)
		catch
			UNDEF
		end
	end

	expression = quote
		import Base: ==, hash
		function =={T<:$typ}(a::T, b::T)
			a === b && return true
			for field in $fields
				fa = $get(a, field)
				fb = $get(b, field)
				# println("fa ", fa, " fb ", fb)
				if !isequal(fa, fb)
					return false
				end
			end
			return true
		end
		function hash(x::$typ, h::UInt)
			for field in $fields
				h = hash($get(x, field), h)
			end
			h
		end
	end
	eval(current_module(), expression)
	nothing
end

end # module

module TestGeneric
	import Generic

	abstract Ta
	type Tb <: Ta
		field1::Integer
		fieldb::String
	end
	type Tc <: Ta
		field1::Integer
		fieldb::String
	end
	immutable Td <: Ta
		field1::Integer
		fieldb::Integer
	end

	Generic.make_equal(Tb) # objects both of type Tb must have all field equal
	Generic.make_equal(Tc) # objects both of type Tb must have all field equal

	ta = Tb(42, "z1")
	te = Tb(42, "z1")
	tb = Tb(42, "z2")
	tc = Tc(42, "z3")
	td = Td(42, 11)

	tlist = [ta, tb, te, tc, td]

	using Base.Test
	@testset "equals wo field1(Te)" begin	
				@testset "ta" begin
			@test ta == ta
			@test ta == te
			@test ta != tb
			@test ta != tc
			@test ta != td
		end
		@testset "te" begin
			@test te == ta && hash(te) == hash(ta)
			@test te == te
			@test te != tb
			@test te != tc
			@test te != td
		end
		@testset "tb" begin
			@test tb != ta
			@test tb != te
			@test tb == tb
			@test tb != tc
			@test tb != td
		end
		@testset "tc" begin
			@test tc != ta
			@test tc != te
			@test tc != tb
			@test tc == tc
			@test tc != td
		end
		@testset "td" begin
			@test td != ta
			@test td != te
			@test td != tb
			@test td != tc
			@test td == td
		end
    end

end # module

using TestGeneric

for T in [:Ta, :Tb, :Tc, :Td]
	eval( :($T = TestGeneric.$T) )
end
for t in [:ta, :tb, :tc, :td, :te]
	eval( :($t = TestGeneric.$t) )
end
