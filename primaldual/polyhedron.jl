# Generate special types of (random) Polyhedrons.

module Polyhedron

# Polyhedron in n dimensions with sphere touching all n*(n+1)/2 1-dimensional edges.
# Edge touching sphere of a Polyhedron.
#
# Necessary and sufficient condition for this property of a polyhedron P[0:n]
# is the extistence of real values x[0:n] such that x[i] + x[j] = norm(P[i] - P[j])
# for all i < j in 0:n.
#
#
# The result ist stored in a Matrix P[1:n,0:n] meaning P[i] = P[:,i] for i = 0:n
# Start with standard coordinates, P[0] = 0, P[i][k] = 0 for all k > i and i in 1:n
# P[:,0] is not stored at all.
# P is an upper triangular matrix.
#
function touching_edges{T}(P::AbstractArray{T,2}, fx::Function)
  nn = size(P,1)
  z = zero(T)
  x0 = fx(0, (z, z, z))
  x = Array{T}(nn)
  for n = 1:nn
    # Assume P[:,k] is already determined for k < n; find P[:,n].
    u = Array{T}(n-1)
    v = Array{T}(n-1)
    for k = 1:n-1
      pkk = P[k,k]
      pk = view(P, 1:k-1, k)
      u[k] = (pkk + (-(x[k]^2-x0^2) + vecdot(pk, pk) - vecdot(pk, u[1:k-1]) * 2) / pkk) / 2
      v[k] = -((x[k]-x0) + vecdot(pk, v[1:k-1])) / pkk
    end
    alpha = x0 ^ 2 - vecdot(u, u)
    beta  = x0     - vecdot(u, v)
	gamma = 1.0    - vecdot(v, v)
    xn = fx(n, (alpha, beta, gamma))
    x[n] = xn
    P[1:n-1,n] = u + v * xn
    P[n,n] = sqrt(alpha + (2 * beta + gamma * xn) * xn)
  end
  P, [x0; x]
end

# Given f(x) = alpha + 2 *beta * x + gamma * x^2
# determine randomly an x > 0 with f(x) > 0 if possible
#
function fx{T}(n::Int, tup::Tuple{T,T,T})
  z = zero(T)
  xmin = z
  xmax = T(Inf)
  
  alpha, beta, gamma = tup
  # println("n = $n, alpha = $alpha, beta = $beta, gamma = $gamma :")
  if gamma == z
    if beta > z
      xmin = max(z, -alpha / beta / 2)
    elseif beta == z
      if alpha < z
        xmin = xmax
      end
    else # beta < z
      if alpha > z
        xmin = z
        xmax = max(-alpha / beta / 2, z)
      else # alpha <= z
        xmin = xmax
      end
    end
  else # gamma != z
    x1 = - beta / gamma
    f1 = alpha + x1 * beta
    if gamma > z
      if f1 < z
        sq = sqrt(x1 ^ 2 - alpha / gamma)
        xmin = x1 + sq
      end
    else # gamma < z
      if f1 < z
        xmin = xmax
      else # f1 >= z
        sq = sqrt(x1 ^ 2 - alpha / gamma)
        xmin = max(z, x1 - sq)
        xmax = max(z, x1 + sq)
      end
    end
  end
  between(xmin, xmax)
end

using Distributions

beta = Beta(2, 2)
gamma = Gamma(2)
# Random value in the given interval. Upper limit may be infinite.
#
function between{T<:Real}(xmin::T, xmax::T)
    mue = T(0.5) # mean value of distribution
	if xmax < T(Inf)
      x = xmin + rand(beta) * (xmax - xmin)
    else
      x = xmin + mue * rand(gamma) / 2
    end
	# x = mue
	println("x = $x ϵ [$xmin, $xmax]")
	x
end

# Given polyhedron by points P[:,1] ... P[:,n]
# calculate symmetric matrix of distances P[i]-P[j]
#
function distances{T}(P::AbstractArray{T,2})
  n, m = size(P)
  m == n + 1 || error("need P[1:n,1:n+1]")

  D = zeros(T, m, m)
  for k = 1:m
    for i = 1:k-1
      D[i,k] = norm(P[:,i] - P[:,k])
    end
  end
  D + D'
end

# Given polyhedron by points P[:,1] ... P[:,n]
# calculate tangent distances x[1:n]
# with P[i] + x[i] * (P[j]-P[i]) / norm(P[i]-P[j]) on edge i-j is tangent to
# edge touching sphere, if taht exists.
#
function tangent_distances{T}(P::AbstractArray{T,2})
  DIK = distances(P)
  n, m = size(P)
  DK  = vec(sum(DIK, 2)) / n
  D = sum(DK) / m
  X = (DK - D) * n / (n-1) + D / 2
  X
end #tangent_distances

end # module




































