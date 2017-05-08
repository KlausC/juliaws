function DownSum{S}(p::AbstractArray{S,})
#DOWNSUM      Rounded downwards result of sum(p)
#
#   resD = DownSum(p)
#
#On return, resD is sum(p) rounded downwards, also in the presence
#  of underflow. Input vector p may be single or double precision.
#
#Implements Algorithm 7.5 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear. 
#Requires (4m+4k+4)n flops for m,k executions of the repeat-until loops
#  in the calls of TransformK, respectively, and not huge dimension.
#
# written  03/03/07     S.M. Rump
#

  nmax = 2^div(Base.significand_bits(S)+1, 2) - 2

  if length(p) <= nmax                      # not huge dimension
    res,R,p,sigma,Ms = TransformK(p,0)      # s-res = R+sum(p)
    delta = TransformK(p,R,sigma,Ms)        # delta faithful rounding of s-res
  else                                      # huge dimension
    res,tau1,tau2,tau,p = AccSumHugeN(p)    # s=tau1+tau2+sum(tau)+sum(p)
    delta = AccSumHugeN([tau1;tau2;tau;p;-res])
  end
  if delta < 0                              # s < res
    res = prevfloat(res)
  end
  res
end
