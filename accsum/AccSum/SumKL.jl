function SumKL{S<:Real}(p::AbstractArray{S,}, K, L = 1)
#SUMKL        Summation 'as if' computed in K-fold precision and stored in L results
#
#   res = SumKL(p,K,L)
#
#On return, sum(res) approximates sum(p) with accuracy as if computed 
#  in K-fold precision, where res comprises of L elements. 
#  Default for L is 1.
#
#Implements algorithm SumKL from
#  S.M. Rump: Inversion of extremely ill-conditioned matrices in floating-point,
#    submitted for publication in JJIAM, 2008.
#
#Reference implementation! Slow due to interpretation!
#

# written  06/23/08     S.M. Rump
#

  n = length(p)
  for i = 1:K-L
    p = VecSum(p)
  end
  res = zeros(S,1,L)
  for k = 0:L-2
    p[1:n-k] = VecSum(p[1:n-k])
    res[k+1] = p[n-k]
  end
  res[L] = sum(p[1:n-L+1])
  res
end
