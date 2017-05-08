#SINKHORNKNOPP normalises a square matrix to be doubly stochastic
#   M = SINKHORNKNOPP(A) takes a nonnegative NxN matrix A and normalises it
#   so that the sum of each row and the sum of each column is unity. M is
#   equal to DAE where D and E are diagonal matrices with positive values
#   on the diagonal.
#
#   M = SINKHORNKNOPP(A, NAME1, VALUE1, ...) allows other parameters to be
#   set. These are:
#
#       'tol' - a positive scalar, default sqrt(EPS)*N. The value is the
#       maximum error in the column sums of M; the row sums will be correct
#       to within rounding error.
#
#       'maxiter' - a positive integer, default 1M. The maximum number of
#       iterations to carry out.
#
#   [M, R, C] = SINKHORNKNOPP(...) also returns normalising vectors R and C
#   such that M = diag(R) * A * diag(C) to within a small tolerance.
#
#   Example
#   -------
#
#       a = toeplitz(1:6);
#       m = sinkhornKnopp(a);
#       disp('a'); disp(a);
#       disp('m'); disp(m);
#       disp('Row and column sums'); disp(sum(m,1)); disp(sum(m,2));
#
#   Convergence
#   -----------
#
#   The algorithm will converge for positive matrices, but may not converge
#   if there are too many zeros in A, depending on their distribution. In
#   such cases it may be necessary to set 'MaxIter' and to check the column
#   sums of M. For some applications adding a small constant to A is
#   recommended.
#
#   Algorithm
#   ---------
#
#   The Sinkhorn-Knopp algorithm, also known as the RAS method and
#   Bregman's balancing method, is used. The code is modified from Knight
#   (2008), avoiding the matrix transpose.
#
#   Reference
#   ---------
#
#   Philip A. Knight (2008) The SinkhornKnopp Algorithm: Convergence and
#   Applications. SIAM Journal on Matrix Analysis and Applications 30(1),
#   261-275. doi: 10.1137/060659624
#   Copyright David Young 2015
#   Translated to Language julia Klaus Crusius (2016)

function sinkhornKnopp{T<:Number}(A::AbstractArray{T,2}; tol::Float64 = 0.0, maxiter = 1000000)
# returning A, r, c
# Input parameter parsing and checking

N, M = size(A)
if N != M
  error("Matrix is not square")
end

if tol == 0.0
    tol = sqrt(eps()) * N
end

if tol <= 0.0
  error("tolerance must be positive")
end

# first iteration - no test
iter = 1
c = 1 ./ sum(abs(A), 1)
r = 1 ./ (A * c.')

# subsequent iterations include test
while iter < maxiter
    iter = iter + 1
    cinv = r.' * A
    # test whether the tolerance was achieved on the last iteration
    if  norm(abs(cinv .* c - 1), Inf) <= tol
        break
    end
    c = 1 ./ cinv
    r = 1 ./ (A * c.')
end

A = A .* (r * c)
A, r, c
end

