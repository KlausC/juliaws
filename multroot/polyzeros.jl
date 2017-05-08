
# PolyZeros represents a set of distinct (potential) roots z[j] of a polynomial
# The z-values are ordered by absolute value to obtain better numeric stability
# when calculating products of the form  Prod over j (x-z[j]) ^mult[j]
# where l defines a pejorative multiplicity structure of the polynomial.
# The stored values are obtaind by z[jj] and mult[jj] from the original arrays.
#
if ! isdefined(:PolyZeros)
immutable PolyZeros{S<:Number}
    z::AbstractVector{S}		# the distinct root values, already ordered
    perm::AbstractVector{Int}	# the permutation z = zorig[perm] maps original ordering
    mult::AbstractArray{Int}	# the multiplicity structure of the zeros, reordered
    d::Dict{Tuple{Int64, Array{Int64,1}}, Array{S,1}}	# memorize intermediate results
end
end

# Constructor for PolyZeros
# Perform initial ordering and set up dictionary
function PolyZeros{S}(zz::AbstractVector{S}, ll::AbstractVector{Int} = Int[])
    m = length(zz)
	if length(ll) == 0 ll = ones(Int,m) end
    if length(ll) != m
        throw(ArgumentError("roots and multiplicities have differing sizes $m and $(length(ll))"))
    end
    jj = sortperm(abs(zz), lt=lessthan)
    d = Dict{Tuple{Int64,Array{Int64,1}}, Array{S,1}}()
    PolyZeros(zz[jj], jj, ll[jj], d)
end

function lessthan{S<:Number, T<:Number}(a::S, b::T)
    ra, rb = real(a), real(b)
    ia, ib = imag(a), imag(b)
    ra < rb ? true : ra > rb ? false : abs(ia) < abs(ib) ? true : abs(ia) > abs(ib) ? false: ia < ib
end

# Determine type of polynomial coefficients given the roots (with type)
function coeffstype{S}(z::AbstractVector{S})
    if isreal(z) return S end
    zc = z[find(x-> !isreal(x), z)]
    zd = zc
    for c = zc
        jj = find(x-> x == conj(c), zd)
        if length(jj) >= 1
            j = jj[1]
            zd = zd[[1:j-1;j+1:end]]
		else
            return S
        end
    end
	return real(S)
end

# Calculate product over j of (x-z[j])
function productmono{S}(z::AbstractVector{S})
    n = length(z)
    c = zeros(S, n+1)
    c[1] = 1
    for j = 1:n
        zj = z[j]
        for i = j:-1:1
            c[i+1] -= zj*c[i]
        end
    end
    isreal(c) ? real(c) : c
end

#
# calculate product over j of (x - z[j])^ll[j]
#
function products{S}(z::PolyZeros{S}, ll::AbstractVector{Int})
    # the greatest potency of 2 which is <= argument
    lastpot(x) = 2^searchsortedlast(2.^(1:62), min(x-1,2^62)+1)
	
    n = sum(ll)	
    s = [1; zeros(S, n)]
	n = 1
	l = ll[z.perm]
    lmax = maximum(l)
	while lmax > 0
        lp = lastpot(lmax)	# greatest potency of 2 <= lmax
		jj = find(k-> k >= lp, l)
		t = get!(z.d, (1,jj)) do
                productmono(z.z[jj])
        end
        k = 2
        while k <= lp
            t = get!(z.d, (k,jj)) do
                    conv(t, t)
            end
			k *= 2
		end
        conv!(s, n, t)
        n += length(t) - 1
		l[jj] -= lp
        lmax = maximum(l)
    end
    isreal(s) ? real(s) : s
end
products{S}(z::PolyZeros{S}) = products(z, z.mult[invperm(z.perm)])

#
# calculate product over j of (x - z[j])^ll[j]
#
function evaluate{S,T}(z::PolyZeros{S}, x::T, ll::AbstractVector{Int})
	
	l = ll[z.perm]
	s = one(promote_type(T,S))
	n = length(ll)
	i = 1
	zz = z.z
	while i <= n
	    zi = x - zz[i]
		li = l[i]
	    if i < n && zz[i+1] == conj(zz[i])
		    m = min(li, l[i+1])
			if m > 0
		        s = s * (real(zi)^2 + imag(zi)^2) ^ m
			end
			if li > m
			    s = s * zi ^ (li - m)
			else
			    s = s * conj(zi) ^ (l[i+1] - m)
			end
			i += 2
		else
	        s = s * zi ^ li
			i += 1
		end
    end
    isreal(s) ? real(s) : s
end
evaluate{S,T}(z::PolyZeros{S}, x::T) = evaluate(z, x, z.mult[invperm(z.perm)])

function invperm(p::AbstractVector{Int})
   n = length(p)
   r = zeros(Int, n)
   r[p] = 1:n
   r
end
