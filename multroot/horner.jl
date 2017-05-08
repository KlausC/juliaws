#

#
# Horner initially stores a list of coefficients of a polynomial, a evaluation point,
# and optionally a corresponding list of absolute errors of the coefficients.
# On demand it modifies the values in order to obtain the values of the selected derivatives
# at the given point and the corresponding error value.
#

if ! isdefined(:Horner)
type Horner{S<:Number,T<:Number,U<:AbstractFloat,V<:AbstractFloat}
    x::S							# evaluation spot
    coeff::AbstractVector{T}		# polynomial coefficients
    devx::U							# absolute error estimatin for x
    devcoeff::AbstractVector{V}		# absolute error estimatins for coefficients
    valid::Int						# evaluation status
end
end

# Create Horner table for evaluation at point x.
# Provide polynomial coefficients and optionally absolute errors for x and the coefficients. 
function Horner{S<:Number,T<:Number,U<:AbstractFloat,V<:AbstractFloat}(x::S, coeff::AbstractVector{T}, devx::U, devcoeff::AbstractVector{V})

    if devcoeff == []
        devcoeff = abs(coeff) * eps(promote_type(T,Float64))
    else
        devcoeff = copy(devcoeff)
    end
    if devx < 0
        devx = abs(x) * eps(real(S))
    end
    if length(coeff) != length(devcoeff) throw(ArgumentError("devcoeff length wrong")) end
    R = promote_type(S,T)
	Horner(x, map(R,coeff), devx, devcoeff, -1)
end

# Evaluate Horner table up to and including m deviations.
# return the m'th deviation of the polynomial and the error estimation
# If m >= degree, all deviations ot the polynomial have been calculated at x.
# Further calls just return the cached results.
function call{S,T,U,V}(h::Horner{S,T,U,V}, m::Int)
#function (h::Horner{S,T,U,V})(m::Int)
    if m < 0 throw(ArgumentError("index $m is not >= 0")) end
	n = length(h.coeff) - 1
    if m > h.valid
        if h.x != 0
            ep = eps(real(T))
            b, e = h.coeff, h.devcoeff
            x, tol = h.x, h.devx
            ax = abs(x)
            for k = max(h.valid+1,0):min(m, n-1)
                for j = 1:n-k
                    bjp = b[j+1] + b[j] * x
                    b[j+1] = bjp
                    e[j+1] = hypot(hypot(abs(e[j+1]), abs(b[j])*tol +ax*e[j]), abs(bjp)*ep)
                end
            end
        end
        h.valid = m
    end
    m <= n ? (h.coeff[n+1-m], h.devcoeff[n+1-m]) : (zero(T), zero(V))
end

#
# evaluate polynomial f and all derivatives at given approximate root-values z
#
function horner_analysis{S<:Number,T<:Number}(z::AbstractVector{S}, f::AbstractVector{T})
    thresh = 1000.0
    m = length(z)
    n = length(z)-1
    if m > length(f)-1
        throw(ArgumentError("polynomial of degree $n with $m zeros"))
    end
	h = Array{Horner}(m)
    for i = 1:m
        h[i] = Horner(z[i], f, eps(abs(z[i])) * 100, map(eps, f) * 100)
    end
	for i = 1:m
        println("root $i: $(z[i])")
        mult = 0
        search = true
		for k = 0:n
            hv, he = h[i](k)
            deviation = abs(hv / (he + eps(zero(real(S)))))
            @printf("%d %d %6.2g %8.8g %8.8g\n", i, k, deviation, abs(hv), he)
            if search && deviation < thresh
                mult = k + 1
            else
                search = false
            end
        end
        println("root $i: $(z[i]) estimated multiplicity: $mult")
    end
end


