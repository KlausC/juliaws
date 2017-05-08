function UpSum{S}(p::AbstractArray{S,})
#UPSUM        Rounded upwards result of sum(p)
#
#   resD = DownSum(p)
#
#On return, resD is sum(p) rounded upwards, also in the presence
#  of underflow. Input vector p may be single or double precision.
#
#Adapted from Algorithm 7.5 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear.
#Requires (4m+4k+4)n flops for m,k executions of the repeat-until loops
#  in the calls of TransformK, respectively.
#
# written  03/03/07     S.M. Rump
#

  nmax = 2^div(Base.significand_bits(S)+1,2) - 2

  if length(p) <= nmax                      # not huge dimension
    res,R,p,sigma,Ms = TransformK(p,0)      # s-res = R+sum(p)
    delta = TransformK(p,R,sigma,Ms)        # delta faithful rounding of s-res
  else                                      # huge dimension
    res,tau1,tau2,tau,p = AccSumHugeN(p)    # s=tau1+tau2+sum(tau)+sum(p)
    delta = AccSumHugeN([tau1;tau2;tau;p;-res])
  end
  if delta > 0                              # s > res
    res = nextfloat(res)
  end
  res
end
