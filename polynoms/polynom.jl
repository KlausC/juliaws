
import Polynomials
using Polynomials
import Base.eps

eps{T}(x::Type{Complex{T}}) = eps(T)
eps{T}(x::Complex{T}) = eps(T)

# shift transformation of variable using Horner-scheme:
# pnew(z) = p(z*s+a) - e.g. Taylor development at point a with additional scaling factor s
#
function polytransform{T<:Number,S<:Number,U<:Number}(p::Poly{T}, a::S = 0, s::U = 1)
    R = promote_type(T, S, U)
	b = map(R, p.a)
	n = length(b) - 1
    if a != 0
        for k = 1:n
		    for j = n:-1:k
			    b[j] += b[j+1] * a
		    end
	    end
    end
    if s != 1
        x = one(T)
        for k = 1:n
            x *= s
            b[k+1] *= x
        end
    end
	Poly(b)
end

# shift transformation of variable using Horner-scheme:
# additionally standard deviation of error is calculated
# pnew(z) = p(z*s+a) - e.g. Taylor development at point a with additional scaling factor s
#
function polytransform1{T<:Number,S<:Number,U<:Number}(p::Poly{T}, a::S = 0, s::U = 1; tol=0)
    R = promote_type(T, S, U)
	b = map(R, p.a)
	e = abs(b) * tol    
	n = length(b) - 1
    if a != 0
        for k = 1:n
		    for j = n:-1:k
			    b[j] += b[j+1] * a
				e[j] = hypot(e[j], abs(b[j+1])*tol+abs(a)*e[j+1])
		    end
	    end
    end
    if s != 1
        x = one(T)
        for k = 1:n
            x *= s
            b[k+1] *= x
			e[k+1] = hypot(e[k+1]*x, b[k+1]*tol*x*n) 
        end
    end
#	println("stddev of result: $e")
	b, e
end

# Interpolate using Newton-method of finite differences.
# given x,f of equal length n+1 create polynomial
# which interpolates y at x[1..n]a
# Alternative to using polyfit(x, y, sym), which is sometimes numerically better
# sometimes worse, but has O(n³) effort, while Newton requires O(n²) operations.
#
function polyinterpol_newton{S, T}(x::AbstractArray{S}, y::AbstractArray{T}, sym::Symbol = :x, R::Type = Float64)
	m::Int64 = length(x)
    m == length(y) || throw(DomainError)
		
	R = promote_type(R, S, T)
	x = Array{R}(x)
	b = Array{R}(y)

	for k::Int64 = 2:m
		for j::Int64 = m:-1:k
			b[j] = (b[j] - b[j-1]) / (x[j] - x[j-k+1])
		end
	end
	n::Int64 = m - 1
	for k::Int64 = 1:n
		for j::Int64 = n:-1:k
			b[j] -= b[j+1] * x[j-k+1]
		end
	end
	Poly(b, sym)
end

# Interpolation using divide and conquer algorithm (own development).
#

function polyinterpol{S, T}(
							x::AbstractArray{S},
							y::AbstractArray{T},
							sym::Symbol = :x,
							R::DataType = Float64)
	m::Int64 = length(x)
	m == length(y) || throw(DomainError)

	R = promote_type(R, S, T)
	xp = Array{R}(x)
	yp = Array{R}(y)
	# println("R == $(R) xp::$(typeof(xp)) yp::$(typeof(xp))")
	polyinterpol_intern(xp, yp, sym)
end

function polyinterpol_intern(x, y, sym::Symbol=:x)
	m::Int64 = length(x)
	
	if m <= 1
		Poly(y, sym)
	elseif m == 2
		a2 = (y[2] - y[1]) / (x[2] - x[1])
		a1 = (y[1]*x[2] - y[2]*x[1]) / (x[2] - x[1])
		Poly( [a1, a2], sym)
	else
		x1, x2 = x[1:2:m], x[2:2:m]
		y1, y2 = y[1:2:m], y[2:2:m]
		p1 = polyinterpol_intern(x1, y1, sym)
		# println("x1 ", x1, "y1 ", y1, "p1", p1)
		q1 = poly(x1, sym) # polynom with roots x1
		# println("q1 ", q1)
		y2 = (y2 - p1(x2)) ./ q1(x2)
		p2 = polyinterpol_intern(x2, y2, sym)
		# println("x2 ", x2, "y2 ", y2, "p2", p2)
		p2 * q1 + p1
	end
end

# Remez - algorithm for Chebyshev approximation
#
function remez(f, df, ddf, n, a=-1.0, b=1.0; tol=1e-30)
	
	xn = sin(linspace(-pi/2, pi/2, n+2)) * ((b-a) / 2) + (a+b)/2
	y = (-1) .^ (0:n+1) * 1.0
	m = n ÷ 2
	yr = [y[1:m-1]; y[m+1:end]]
	x = zeros(xn)

	while true
		x = xn
		xr = [x[1:m-1]; x[m+1:end]]
		xm = x[m]
		p = polyinterpol(xr, f(xr))	
		q = polyinterpol(xr, yr)

		E = (p(xm) - f(xm)) / (q(xm) - y[m])
		pmq = p - q * E

		xn = x - (df(x) - polyder(pmq)(x)) ./ (ddf(x) - polyder(pmq, 2)(x))
		xn[1] = x[1]
		xn[n+2] = x[n+2]

		println("x  ", x')
		println("xn ", xn')

		if all(abs(x-xn) .<= tol) return E, xn, pmq end
	end	

end	

# gcdroot{T}(p::Poly{T}, tol = eps(T)) = gcdroot(reverse(p.a), tol)

