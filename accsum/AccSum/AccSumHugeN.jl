function AccSumHugeN{S}(p::AbstractArray{S,})
# [res,tau1,tau2,tau,p]
#ACCSUMHUGEN  Faithful rounding of sum(p) for huge dimension
#
#   res = AccSumHugeN(p)
#
#On return, res is a faithful rounding of sum(p), also in the presence
#  of underflow. Input vector p may be single or double precision.
#
#Implements Algorithm 8.1 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear. 
#Requires (4m+4K+3)n flops for m executions of repeat-until loop in 
#  Transform and K executions of the while-loop.
#
# written  03/03/07     S.M. Rump
#

  nmax = 2^(Base.significand_bits(S)-2)
  rm = realmin(S)
  if length(p) > nmax
    error("maximum length of input vector for AccSumHugeN $(nmax).")
  end
 
  kPhi = 2
  tau1,tau2,p,sigma,Ms = Transform(p, 0, kPhi)   # Ms = 2^M
  if sigma <= rm                # p_i identical zero
    res = tau1
    tau = zeros(S,1,1)          # make sure tau has same precision
    return res, tau1, tau2, tau, p
  end
  u = 0.5*eps(S)
  tau = zeros(S,1,16)
  K = 0
  phi = Ms*u
  factor = 2 * Ms * Ms * u
  sigmas = phi * sigma

  while true
    K = K+1
    sigma = sigmas
    tau[K],p = ExtractVector(p, sigma)
    sigmas = phi * sigma
    if ( factor * sigma <= abs(tau1) ) || ( sigma <= rm )
      taus = tau2 + sum(p)
      for k = K:-1:1
        taus = taus + tau[k]
      end
      res = tau1 + taus
      return res, tau1, tau2, resize!(tau), p
    end
	if K >= length(tau)
	  resize!(tau, K*2)
	end
  end
end
