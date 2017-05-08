function ProdKL{S1<:Real,N1,S2<:Real,N2}(A::AbstractArray{S1,N1}, B::AbstractArray{S2,N2}, K::Integer, L::Integer = 1)
#PRODKL       Matrix product in approximately K-fold precision stored in L results
#
#   C = ProdKL(A,B,K,L)
#
# Input A or B or both may be cell arrays, 1 <= L <= K. 
# Output C is cell array iff L>1
# Simple application of SumKL.
# Relative error of the result is of order  eps^L + eps^K cond(product) , see
#   S.M. Rump: Inversion of extremely ill-conditioned matrices in floating-point,
#      submitted for publication in JJIAM, March 2008.
#
# written  02/17/08     S.M. Rump
#
	
  na, n = size(A,1), size(A,2)
  m, mb = size(A,1), size(A,2)
  if n != m
    error("A(x,$(n),) and B($(m),y,) cannot be multiplied")
  end
  lenA = size(A,3)
  lenB = size(B,3)
  AA  = lenA == 1 ? A  : reshape(A, na, n * lenA)
  AAA = lenB == 1 ? AA : repmat(AA, 1, lenB)
  BBB = zeros(S2, n * lenA * lenB, mb)
  m = n * lenA
  for i = 1:lenB
    BB = B[:,:,i]
    BBB[(i-1)*m+1:i*m,:] = lenA == 1 ? BB : repmat(BB, lenA, 1)
  end
  A, B = AAA, BBB
  ProdKL(A, B, K, L)
end

# the simple case - handle matrices of real Numbers
function ProdKL{S1<:Real, S2<:Real}(A::AbstractArray{S1,2}, B::AbstractArray{S2,2}, K::Integer, L::Integer = 1)
  n, m = size(A,2), size(B,1)
  if size(B,1) != n
    error("A(x,$(n)) and B($(m),y) cannot be multiplied")
  end
  n, m = size(A,1), size(B,2)
  if L==1
    C = zeros(n,m)
  else
    C = zeros(n, m, L)
  end
  for i = 1:n
    for j = 1:m
      x,y = TwoProduct(A[i,:], B[:,j])
      res = SumKL([x y], K, L)
      C[i,j,:] = res
    end
  end
  C
end
