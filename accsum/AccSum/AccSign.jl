function AccSign{S}(p::AbstractArray{S,})
#ACCSIGN      Computes sign(sum(p))
#
#   S = AccSign(p)
#
#On return, S is the sign of sum(p), also in the presence
#  of underflow. Input vector p may be single or double precision.
#
#Implements Algorithm 8.2 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear. 
#Requires (4m+2)n flops for m executions of repeat-until loop in Transform.
#
# written  03/03/07     S.M. Rump
#

  nmax = 2^Base.significand_bits(S)
  if length(p) > nmax
    error("maximum length of input vector for AccSign $(int2str(nmax)).")
  end

  kPhi = 1
  tau1,tau2,p = Transform(p,0,kPhi)
  sign(tau1)
end
