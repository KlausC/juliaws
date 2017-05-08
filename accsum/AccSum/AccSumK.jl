function AccSumK{S<:Real}(p::AbstractArray{S,}, K::Integer)
#ACCSUMK      Res represents K-fold faithful rounding of sum(p)
#
#   Res = AccSumK(p,K)
#
#On return, Res represents a K-fold faithful rounding of sum(p), also in the 
#  presence of underflow. Input vector p may be single or double precision.
#
#Implements Algorithm 6.4 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear. 
#Requires (4m+5K+3)n flops for m executions of repeat-until loop in the
#  first and one execution of the repeat-until loops in subsequent calls
#  of TransformK.
#
# written  03/03/07     S.M. Rump
#

  nmax = 2^div(Base.significand_bits(S)+1,2) - 2
  rm = realmin(S)
  if length(p) > nmax
    error("maximum length of input vector for AccSumK $(int2str(nmax)).")
  end

  R = 0
  Res = zeros(S, 1, K)
  Res[1],R,p,sigma,Ms = TransformK(p, R)
  if abs(Res[1]) <= rm
    return Res
  end
  for k = 2:K
    Res[k],R,p,sigma,Ms = TransformK(p, R, sigma, Ms)
    if abs(Res[k]) <= rm
      return Res
    end
  end
  Res
end
