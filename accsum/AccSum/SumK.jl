function SumK(p::AbstractArray, K)
#SUMK         Summation 'as if' computed in K-fold precision
#
#   res = SumK(p,K)
#
#On return, res approximates sum(p) with accuracy as if computed 
#  in K-fold precision.
#
#Implements algorithm SumK from
#  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
#    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
#Requires (6K-5)n flops.
#
# written  03/03/07     S.M. Rump
#
  n = length(p)
  for k = 1:K-1
    for i = 2:n
      p[i],p[i-1] = TwoSum(p[i],p[i-1])
    end
  end
  sum(pi[1:end-1]) + p[end]
end
