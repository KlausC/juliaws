
# calculate the probabiliy, that in a set of n persons all 
# have birthday at different dates of the year.
#

# More abstract the problem is as follows:
# There is a set of n slots and a number k of undistinguishable tokens.
# For each slot i, there is a probability p[i] for a token to be assigned to it.
# After all tokens are assigned, the number a[i] of tokens assigned to slot i is known.
# Sorting the array a and filtering out zeros defines a pattern a[:]. sum a[i] = k. 
# Calculate the probability prob(a, p) for each pattern for n slots.
# 

module birthday

# simple case: uniform distribution p[i] = 1/n for i = 1:n, a[j] = 1 for j = 1:k.
function prob(k::Int, n::Int = 365)
  s = 1.0
  for j = 1:k-1
    s *= (n - j) / n
  end
  s
end

# remember all function evaluations
#
function memo(f::Function)
  vals = Dict()
  function g(args::Any...)
    get!(vals, args) do
      ff = f(args...)
      # println("f($args) = $ff")
      ff
    end
  end
  g
end

# momentum function for a given discrete distribution
#
function momentum{T<:Real}(p::Vector{T})
  m1 = sum(p)
  p = p ./ m1
  mom = Dict{Int,eltype(p)}()
  function m(k::Int)
    get!(mom, k) do
      k == 1 ? one(eltype(p)) : sum(p .^ k)
    end
  end
  m
end

# probability of k persons have distinct birthdays when
# birthdays are distributed according to p.
#
function prob{T}(p::Vector{T}, k::Int)
  vm = getvm(p)
  a = ones(Int, k)
  vm(a)
end

# generate a function pm(a) giving the probability for a pattern
# p(a) is independent of the sorting of a.
#
function getvm{T}(p::Vector{T})
  mom = momentum(p)
  n = length(p)
  function v(a::Vector{Int})
    l = length(a)
    l == 0 && return one(eltype(p))
    l == 1 && return mom(a[1])
    al = a[l]
    b = a[1:end-1]
    s = mom(al) * vm(b)
    for i = 1:l-1
      c = copy(b)
      c[i] += al
      s -= vm(sort(c, rev=true))
    end
    s
  end
  vm = memo(v)
  function pm(a::Vector{Int})
    (n - length(a) + 1 <= 0) && return 0
    k = sum(a)
    p = countpattern(n, a)
    for i = n:-1:n-length(a)+1
      p ÷= i
    end
    p * vm(a)
  end
  pm
end

#special case p[i] = 1 / n for i = 1:n
#
function getvm(n::Integer)
  pp = 1 // BigInt(n)
  mom(k::Int) = pp ^ (k-1)
  function v(a::Vector{Int})
    l = length(a)
    l == 0 && return one(eltype(p))
    l == 1 && return mom(a[1])
    s = mom(sum(a))
    for k = 2:l
      s *= (n-k+1)
    end
    s
  end
  vm = memo(v)
  function pm(a::Vector{Int})
    (n - length(a) + 1 <= 0) && return 0
    k = sum(a)
    p = countpattern(n, a)
    for i = n:-1:n-length(a)+1
      p ÷= i
    end
    p * vm(a)
  end
  pm
end

# Determine all summations of positive integers giving n.
# The greatest addend is k.
#
function allsums(n::Int, k::Int)
  all = Vector{Vector{Int}}()
  n <= 0 && return all
  k >= n && append!(all, [[n]])
  for i in min(k,n):-1:1
    b = allsums(n-i, i)
    append!(all, map(bb->[i; bb], b))
  end
  all
end
    
# Determine count of all summations of positive integers giving n.
# The greatest addend is k.
#
function countallsums(n::Int, k::Int)
  all = Int128(0)
  n <= 0 && return all
  (k >= n) && (all += 1)
  for i in min(k,n):-1:1
    all += countallsums(n-i, i)
  end
  all
end

# count patterns
# Generate all k-tuples of values 1:n
# Count all patterns of muliplicites of entries.
# for example (3,2,3,1) counts for pattern (2,1,1), because
# 3 appears twice, 1 and 2 each once.
function countpatterns_obsolete(n::Int, karg::Int)
  function next!(a::Vector{Int}, n::Int)
    a[end] += 1
	for k = length(a):-1:2
	  if a[k] >= n
	    a[k] = 0
		a[k-1] += 1
	  else
	    break
	  end
    end
  end
  a = zeros(Int, karg)
  pd = Dict{Vector{Int}, Int}()
  while a[1] < n
    d = Dict{Int,Int}()
    for i = 1:karg
	  ai = a[i]
	  d[ai] = count(aa->(ai == aa), a)
    end
	s = sort(map(identity, values(d)), rev=true)
    pd[s] = get(pd, s, 0) + 1
    next!(a, n)
  end
  pd
end

function countpatterns(n::Int, k::Int, k2 = k)
  Dict(map(aa->(aa, countpattern(n, aa)), allsums(k, k2)))
end

function countpattern(n::Int, a::Vector{Int})
  k = sum(a)
  k == 0 && return BigInt(1)
  a1 = a[1]
  BigInt(n) * binomial(BigInt(k), a1) * countpattern(n-1, a[2:end]) ÷ count(x->(x==a1), a) 
end

function countpattern(n::Int, a::Vector{Int}, cum::Rational{BigInt})
  k = sum(a)
  k == 0 && return cum
  a1 = a[1]
  countpattern(n-1, a[2:end], cum * n * binomial(k, a1) // count(x->(x==a1), a)) 
end

end # module

bd = birthday
