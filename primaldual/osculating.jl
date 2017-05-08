module SimplexCircles

export outer_circle, inner_circle

# Given m+1 points in Rn with m <= n forming a m-dimensional simplex.
# Shift representation with first point as origin.
# Perform orthogonal transformation to work in m-dimensional subspace.
# Return n, m, shifting vector, orthogonal transformation matrix Q in Rn,m (thin form),
# transformed points R in Rm,m.
#
function prepare{T <: Real}(A::Array{T,2})
  n, m = size(A); m -= 1
  0 <= m <= n || error("m $m must be <= n $n and >= 0")

  # shift by one of the column vectors contained in A
  ac = A[:,1]
  B = A[:,2:m+1] - ac * ones(m)'

  # Orthogonal transformation to triangular form
  Q, R = qr(B)
  n, m, ac, Q, R
end

# Given m+1 points in Rn with (m <= n). Find the m-dimensional hypershere in the
# m-dimensional affine hyperspace of Rm, which contains all m+1 points.
# return center and radius of the sphere.
# If the points are denoted Pk = A[:,k] for k = 1:m+1
# then P1 = ac, Pk = ac + Q * R[:,k-1] for k = 2:m+1
#
function outer_circle{T <: Real}(A::Array{T,2})
  n, m, ac, Q, R = prepare(A)
  z = similar(Q, m)
  rho2 = 0.0
  for k = 1:m
    z[k] = ( ( sum((z[1:k-1] - R[1:k-1,k]).^2) - rho2 ) / R[k,k] + R[k,k] ) / 2.0 
    rho2 = z[k] ^ 2 + rho2
  end
  c = Q * z + ac
  # sqrt(rho2) and c are first estimations.
  # try to reduce rounding errors by executing single Newton-iteration. 
  b = (sum((A - c * ones(m+1)') .^ 2, 1))'
  dc = Q * (R' \ (b[2:end] - b[1]))
  dr2 = sum(b - A' * dc) / float(m+1)
  dr2 += vecdot(c, dc)
  sqrt(dr2), c + dc / 2.0
end

# Given m+1 points in Rn with (m <= n). Find the m-dimensional volume in the
# m-dimensional convex hull of the points.
#
function volume{T <: Real}(A::Array{T,2})
  n, m, ac, Q, R = prepare(A)
  vol = prod(diag(R)) / float(m)
  vol
end

# Given m+1 points in Rn with (m <= n). Find the m-dimensional volume and the
# (m-1)-dimensional surface of the m-dimensional convex hull of the points.
#
function inner_circle{T <: Real}(A::Array{T,2})
  n, m, ac, Q, R = prepare(A)
  N = inv(UpperTriangular(R))'
  n0 = N * ones(m)
  D = sqrt(sum(N .^2, 1))
  d0 = norm(n0)
  d = (sum(D) + d0) .\ 1.0
  c = R * D' * d
  d, vec(Q * c  + ac), Q * [ N ./ D * d -n0 / d0 *d ]
end

# given m+1 tangent points and normals of the m+1 facets of a simplex, calculate the
# extreme points of the simplex.
#
function outer_points{T <: Real, U <: Real}(P::Array{T,2}, N::Array{U,2})
  n, m, ac, Q, R = prepare(P)
  NQ = Q' * N
  X = zeros(n,m+1)
  b = zeros(m+1)
  b[2:m+1] = map( k -> vecdot(NQ[:,k+1], R[:,k]), 1:m) 
  X[:,1] = NQ[:,1:m]' \ b[1:m]
  for k = 1:m
    X[:,k+1] = [NQ[:,1:k-1] NQ[:,k+1:m+1]]' \ [b[1:k-1]; b[k+1:m+1]]
  end
  ac .+ Q * X 
end

using Base.Test

# generate m points uniformly distributed on surface of n-dimensional sphere in R^n.
function test_matrix(n::Int, m::Int)
  A = randn(n, m)
  sqrt(sum(A .^ 2, 1)) .\ A
end

# 

# generate m points of R^n with inner circle radius 1.
function test_matrix_inner(n::Int, m::Int)
  P = test_matrix(n, m)
  d, c, N = inner_circle(P)
  -N / norm(N[:,1])
end

function test_inner_circle(n::Int)
  m = n + 1
  N1 = test_matrix_inner(n, m)
  A = outer_points(N1, -N1)
  d, c, N = inner_circle(A)
  d, c
end

@testset "osculating" begin

@testset "outer circle $n" for n in 2:50
function test_outer_circle(n::Int)
  m = n + 1
  A = test_matrix(n, m)
  d, c = outer_circle(A)
end

  d, c = test_outer_circle(n)
  @test isapprox(d, 1.0, rtol = eps()*n*10)
  @test isapprox(norm(c, Inf), 0.0, atol = eps()*n*100)

end # test outer circle

@testset "inner circle $n" for n in 2:50
  d, c = test_inner_circle(n)
  @test isapprox(d, 1.0, rtol = eps()*n*100)
  @test isapprox(norm(c, Inf), 0.0, atol = eps()*n*1000)
end # test inner circle
end # test osculating
end # module
