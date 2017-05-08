
using ValidatedNumerics
#
# HornerInt initially stores a list of coefficients of a polynomial, a evaluation point,
# and optionally a corresponding list of absolute errors of the coefficients.
# On demand it modifies the values in order to obtain the values of the selected derivatives
# at the given point and the corresponding error value.
#

if ! isdefined(:HornerInt)
type HornerInt{S<:AbstractFloat, T<:AbstractFloat}
    x::Interval{S}						# evaluation spot
    coeff::AbstractVector{Interval{T}}	# polynomial coefficients
    valid::Int						# evaluation status
end
end

# Create Horner table for evaluation at point x.
# Provide polynomial coefficients and optionally absolute errors for x and the coefficients. 
function HornerInt{S<:Number,T<:Number,U<:AbstractFloat,V<:AbstractFloat}(xdown::S, xup::T, coeffdown::AbstractVector{U}, coeffup::AbstractVector{V})

    if coeffup == []
        coeffup = copy(coeffdown)
    end
    xdown, xup = min(xdown, xup), max(xdown, xup)
    coeffdown, coeffup = min(coeffdown, coeffup), max(coeffdown, coeffup)
    if length(coeffdown) != length(coeffup) throw(ArgumentError("coeffup length wrong")) end
    R = promote_type(S,T,U,V)
	HornerInt(R(xdown), R(xup), map(R,coeffdown), map(R,coeffup), -1)
end

# Evaluate Horner table up to and including m deviations.
# return the m'th deviation of the polynomial and the error estimation
# If m >= degree, all deviations ot the polynomial have been calculated at x.
# Further calls just return the cached results.
function call{S}(h::HornerInt{S}, m::Int)
    if m < 0 throw(ArgumentError("index $m is not >= 0")) end
	n = length(h.coeffdown) - 1
    if m > h.valid
        if h.xdown != 0 || h.xup != 0
            r = rounding(S)
            setronding(S, RoundDown)
            b, bu = h.coeffdown, h.coeffup
            x, xu = h.xdown, h.xup
            for k = max(h.valid+1,0):min(m, n-1)
                for j = 1:n-k
                    b[j+1], bu[j+1] = intmin((b[j+1], bu[j+1), (b[j], bu[j]), (x, xu))
                end
            end
            setRounding(S, r)
        end
        h.valid = m
    end
    m <= n ? (h.coeffdown[n+1-m], h.coeffup[n+1-m]) : (zero(S), zero(S))
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


