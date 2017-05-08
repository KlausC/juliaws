
module OptUnits

export optunit

typealias MI AbstractArray{Int,2}
typealias VI AbstractVector{Int}

" count all non-zero entries"
xnorm(b::VI) = Int(norm(b, 0))
ynorm = xnorm

norm1(a::VI) = sum(abs(a))
norm0(a::VI) = sum(ifelse(a .== 0, 0, 1))
normx(a::VI) = 0 + any(x -> x != 0, a)

" find k0, k1 so norm1(a*k-b) is minimal for all k in k0:k1"
function heuristicopt(a::VI, b::VI, k::Any)
	ab = filter(ab -> ab[1] != 0, zip(a,b))
	c = sort(unique([[fld(x[2], x[1]) for x in ab]; [cld(x[2], x[1]) for x in ab]]))
	c = filter(cc -> cc != 0, c)
	function ntup(cc::Int)
		r = a * cc - b
		(norm0(r), norm1(r), normx(r), k, cc)
	end	
	map(ntup, c)
end
		



function optunit(M::MI, b::VI)
	n, m = size(M)
	n == size(b,1) || throw(ArgumentError("dimensions of M and b do not match"))
	mb = sum(abs(M),2)
	prim = sum(k ->b[k] != 0 && mb[k] == 0, 1:n)
	bb = [ifelse(mb[k] == 0, 0, b[k]) for k = 1:n]
	opt, z, res = optunit_rec(M, bb, zeros(Int, m), norm0(bb))
	opt + prim, z, b - M * z
end

function optunit_rec(M::MI, b::VI, z::VI, bound::Int)
	n, m = size(M)
	n0 = norm0(z)
	res = b
	bt = n0 + normx(b)
	if bt >= bound
		# println("1: $b $z $bound => $res $z : >= $bt not accepted")
		return bt, z, res
	end
	opt = norm0(z) + norm0(b)
	if n0 >= m
		println("2: $b $z $bound => $res $z - $bt:$opt full")
		return opt, z, res
	end
	heu = Vector{Tuple}()
	for k = 1:m
		if z[k] == 0
			append!(heu, heuristicopt(M[:,k], b, k))
		end
	end
	sort!(heu)
	zopt = z
	for x in heu
		n0, n1, nx, k, cc = x
		zz = copy(z)
		zz[k] = cc
		bb = b - M[:,k] * cc
		opt2, z2, res2 = optunit_rec(M, bb, zz, min(opt, bound))
		if opt2 < min(opt, bound)
			opt = opt2
			zopt = z2
			res = res2
		end
	end
	println("3: $b $z $bound => $res $zopt - $bt:$opt")
	opt, zopt, res
end

end # module


module SIOptUnits

using SIUnits
using OptUnits
using PhysicalConstantsSI

function matrix()
	M = zeros(Int,length(SIUnits.sidims(Meter)),0)
	for y in [names(SIUnits); names(PhysicalConstantsSI)]
		x = eval(y)
		if isa(x, SIUnits.SIUnit)
			col = collect(SIUnits.sidims(x))
			if OptUnits.norm0(col) > 1
				println("$y $col")
				M = [M col]
			end
		end
	end
	M
end


end # module




