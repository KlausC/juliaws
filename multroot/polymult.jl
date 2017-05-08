import Polynomials
using Polynomials

function polymult{T}(z::AbstractVector{T}, l::AbstractVector{Int})
#
# construct monic polynomial from its roots
# use sorting and a special order of processing to increase accuray.
# special implementation of poly()
# returns Real polynomal, if all conjugate roots with equal multiplicities.
#
#   INPUT:
#            z vector of distinct roots of polynomial
#            l vector of corresponding multiplicities
#   OUTPUT:
#            ff Poly representing the product (x - z[k]) ^ l[k]
#
    # println("polymult($(z), $(l)")
    z = PolyZeros(z)
    products(z, l)
end

# Calculate the coefficient operator of the pejorative manifold
# PolyRoots is used to leverage calculation of (x-z[j]) - products
#
function polyG{S}(z::PolyZeros{S})
    products(z)[2:end]
end

# The derivative of G at point z
function polyDG{S}(z::PolyZeros{S})
    m = length(z.z)
    ll = z.mult
	s = products(z, max(ll-1, 0))
    n = sum(ll)
    Df = zeros(S, n, m)
    l = map(k->1, ll)
    for j = 1:m
        l[j] = 0
        Df[1:end,j] = - conv(s, products(z, l)) * ll[j]  
        l[j] = 1
    end
    Df
end

function conv!{S<:Number,R<:Number}(r::AbstractVector{R}, m::Int, q::AbstractVector{S})
# convolution of two vectors in place- represent polynomial multiplication
#
    n, m1 = length(q), length(r)
    mn = m + n - 1
	if m1 < mn
        r = resize!(R, mn)
        r[m+1:end] = zeros(R, mn-m1)
	end
    Z = zero(R)
    for k = mn:-1:1
		j1 = min(k, m)
        j2 = min(k, n)
        s = Z
        for j = (k+1-j1): j2
            s += r[k+1-j] * q[j]
        end
        r[k] = s
    end
	isreal(r) ? real(r) : r
end
function conv{T<:Number,S<:Number}(p::AbstractVector{T}, q::AbstractVector{S})
    m, n = length(p), length(q)
    mn = m + n - 1
    r = [p; zeros(S, mn-m)] 
    conv!(r, m, q)
end
    
import Polynomials.roots
function roots{T}(pc::AbstractVector{T})
# redirect to implementation in Polynomials
	p = Poly(reverse(pc))
	r = roots(p)
    r
end

function conv{T<:Number,S<:Number,N,M}(p::AbstractArray{T,N}, q::AbstractArray{S,M})
# convolution of two vectors - represent polynomial multiplication
#
	conv(vec(p), vec(q))
end

function deriv{T}(p::AbstractVector{T})

    n = length(p) - 1
    map(k -> (p[k] * (n-k+1)), 1:n)
end
