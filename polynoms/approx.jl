
# approximate a given function f(x) by a parametric function g(x, p)
# by minimizing max(|f(xi) - g(xi, p)|) where xi varies over a set of
# values.
# Use derivative-free optimization method.
# Return tuple of status, minimal max-deviation optimal parameter p, 

import BlackBoxOptim
using BlackBoxOptim

function approx(f, g, xs, startvector, range = (-2.0, 2.0))

	diff(x, p) = abs(f(x) - g(x, p))
    function fun(p)
		h(x) = diff(x, p)
		maximum(map(h, xs))
	end
	
	st, bestvector, pp = generic_optimize(fun, startvector, range)
	gp(x) = g(x, bestvector)
	@vectorize_1arg(Number, gp)
	(st, gp, fun, bestvector, pp)

end

function generic_optimize(fun, par, sr)

	println("start: ", par)
	method = :adaptive_de_rand_1_bin_radiuslimited
	pp = bboptimize(fun; SearchRange = sr, NumDimensions = length(par),
    					 Method = method, MaxSteps = 100000, MinDeltaFitnessTolerance = 1e-16 )

	bestvector = pp.best_candidate
	(pp.stop_reason, bestvector, pp) 
	# pp = optimize(fun, par, method = :nelder_mead, grtol = 1e-32, ftol = 1e-33)
end


function polynom(x::Number, par::AbstractArray{Float64,1})
	op(l,r) = r * x + l
	foldr(op, par)
end

function polynom_even(x::Number, par::AbstractArray{Float64,1})
	xs = x * x
	op(l,r) = r * xs + l
    foldr(op, par) * x
end

function polynom_odd(x::Number, par::AbstractArray{Float64,1})
	xs = x * x
	op(l,r) = r * xs + l
    foldr(op, par) * x
end

function ratfun(x::Number, num::AbstractArray{Float64,1}, denom::AbstractArray{Float64,1})
	op(l,r) = r * x + l
	foldr(op, num) / foldr(op, [one(x); denom])
end

function ratfun_even(x::Number, num::AbstractArray{Float64,1}, denom::AbstractArray{Float64,1})
	xs = x * x
	op(l,r) = r * xs + l
	foldr(op, num) / foldr(op, [one(x); denom])
end

function ratfun_odd(x::Number, num::AbstractArray{Float64,1}, denom::AbstractArray{Float64,1})
	xs = x * x
	op(l,r) = r * xs + l
	foldr(op, num) / foldr(op, [one(x); denom]) * x
end

function approxfunction(what = 'a'; degnom = 0)
	if what == 'a'
		if degnom == 0
			fun = polynom
		else
			fun(x, p) = ratfun(x, p[1:end-degnom], p[end-degnom+1:end])
		end
	elseif what == 'e'
		if degnom == 0
			fun = polynom_even
		else
			degnom % 2 != 0 && throw(ArgumentError("degree of nominator must be even"))
		    fun(x, p) = ratfun_even(x, p[1:end-degnom÷2], p[end-degnom÷2+1:end])
		end
	elseif what == 'o'
		if degnom == 0
			fun = polynom_odd
		else
			degnom % 2 != 0 && throw(ArgumentError("degree of nominator must be even"))
		    fun(x, p) = ratfun_odd(x, p[1:end-degnom÷2], p[end-degnom÷2+1:end])
		end
	else
		throw(ArgumentError("what must not be '$(what)'"))
	end
	fun
end
