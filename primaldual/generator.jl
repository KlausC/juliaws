module SubsetIterators

export NoverK

immutable NoverK
  n:: Int
  k:: Int
end

Base.start(it::NoverK) = 0
Base.done(it::NoverK, s::Integer) = s >= length(it)
Base.next(it::NoverK, s::Integer) = (next_impl(it.n, it.k, s), max(s,0) + 1)

function next_impl(n::Int, k::Int, s::Integer) 
  if s <= 0 
    map(identity, 1:k)
  else
    k0 = binomial(BigInt(n-1), k)
    if s < k0
      next_impl(n - 1, k, s)
    else
      if k > 0
        [ next_impl(n - 1, k - 1, s - k0); Int(n) ]
      else
	    error("k < 0")
	  end
    end
  end
end

Base.eltype(::Type{NoverK}) = Array{Int,1}
Base.length(it::NoverK) = binomial(BigInt(it.n), it.k)

end # module
