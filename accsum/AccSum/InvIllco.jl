function InvIllco{S<:Real, N}(A::AbstractArray{S,N}, iter::Any = nothing)
#INVILLCO     Inverse of extremely ill-conditioned matrices
#
#   R = InvIllco(A,iter)
#
#On return, R is a cell array such that sum(R{i}) is an approximate inverse
#  of A and  I-R*A  is convergent. 
#The input matrix A itself may be a cell array. This is convenient to store
#  not exactly representable input data in higher precision. Then R is an
#  approximate inverse of sum(A{i}).
#
#The parameter iter is optional. If specified, an extra iteration is 
#  executed to produce  I-R*A  of order eps.
#
#Example:
#
#  n = 20;              % dimension os matrix
#  C = 1e50;            % anticipated condition number
#  A = randmat(n,C);    % ill-conditioned matrix
#  R = InvIllco(A);     % approximate inverse, stored in cell array R
#  norm(accdot(R,A,-1,eye(n)),'fro')     
#
#Implements algorithm InvIllco from
#  S.M. Rump: Inversion of extremely ill-conditioned matrices in floating-point,
#    submitted for publication in JJIAM, 2008.
#
# written  06/23/08     S.M. Rump
#
  n = size(A,1)
  if size(A,2) != n || N < 2 || N > 3
    error("InvIllco requires square Matrix A")
  end
  R = inverse(A[:,:,1])        # valid for A[:,:,1] == A if A has 2 dimensions
  Eps = eps(S)
  k = 1
  while true
    k = k + 1
    P = ProdKL(R, A, k, 1)
    X = inverse(P)
    R = ProdKL(X, R, k, k)
    if norm(X) * norm(P) < .01 / Eps # ???
      if iter != nothing iter = nothing else break end
    end
  end
  R
end

# If inv(P) throws SingularException, silently disturb P and try again.
function inverse{S<:AbstractFloat}(P::AbstractArray{S,2})
    sz = size(P)
    if sz[1] != sz[2]
        error("matrix P is not quadratic")
    end
    Eps = eps(S)
    X = nothing
    while X == nothing
        X = try
            inv(P)
        catch ex
            if typeof(ex) != LinAlg.SingularException
                rethrow
            end
            P += map(S, randn(sz)) * Eps .* (abs(P) + ones(S,sz) * realmin(S))
        end
    end
    X
end
