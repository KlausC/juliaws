
module CopyModify

export copy_modify
import Base:	copy

@deprecate copy(x::Union{SimpleVector,Module,Symbol,LambdaInfo,GlobalRef,DataType,Union,Task}) identity(x)

# "make a new instance if mutable data type with fields or return identical"
function Base.copy(x::Any)
	ty = typeof(x)
	nf = nfields(ty)
	(isimmutable(x) || isbits(ty) || nf == 0) && return x 
	copy_mutable(x, ty, nf)
end

copy_modify(x::Union{SimpleVector,Array,Module,Symbol,LambdaInfo,GlobalRef,DataType,Union,Task}, ::Dict) = error("copy_modify of $(typeof(x)) not supported")

function copy_modify(x::Any, dic::Dict{Symbol})
	ty = typeof(x)::DataType
	nf = nfields(ty)
	isempty(dic) && return copy(x)
	fields = fieldnames(x)
	isempty(intersect(keys(dic), fields)) && return copy(x)

	y = ccall(:jl_new_struct_uninit, Any, (Any,), ty)::ty
	unmod = true
	for i in 1:nf
		orig = isdefined(x, i) ? getfield(x, i) : ty
		newv = getval(dic, fields[i], orig)::ty.types[i]
		if newv !== ty # use ty as an unique identifier indicating undefined orig
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
	y
end

@inline function getval(dic::Dict, f::Symbol, val)
	if haskey(dic, f)
		v = dic[f]
		isa(v, Function) ? v(val) : v
	else
		val
	end
end

end # CopyModify
