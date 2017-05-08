# count patterns
# Generate all k-tuples of values 1:n
# Count all patterns of muliplicites of entries.
# for example (3,2,3,1) counts for pattern (2,1,1), because
# 3 appears twice, 1 and 2 each once.
function countpatterns(n::Int, k::Int)
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
  a = zeros(Int, k)
  pd = Dict{Vector{Int}, Int}()
  while a[1] < n
    d = Dict{Int,Int}()
    for i = 1:k
	  ai = a[i]
	  d[ai] = count(aa->(ai == aa), a)
    end
	println("$a: $d")
	s = sort(map(identity, values(d)), rev=true)
    pd[s] = get(pd, s, 0) + 1
    next!(a, n)
  end
  pd
end

