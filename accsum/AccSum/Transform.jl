function Transform{S<:Real}(p::AbstractArray{S,}, rho = 0, kPhi = 0, sigma = 0, Ms = 0)
# [tau1,tau2,p,sigma,Ms]
#TRANSFORM    Error-free transformation of rho+sum(p)
#
#   [tau1,tau2,ps] = Transform(p,rho,kPhi)
#
#On return, tau1+tau2+sum(ps_i) = rho+sum(p_i), also in the presence
#  of underflow. Input vector p may single or double precision.
#
#Parameter rho optional, default 0
#
#Parameter kPhi (default 0) specifies stopping criterion:
#  kPhi  0   2^(2M) eps     for AccSum, TransformK, AccSumK, NearSum
#        1   2^M eps        for AccSign
#        2   eps            for AccSumHugeN
#Length of input vector not checked.
#
#Implements Algorithm 4.1 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
#    Faithful Rounding, to appear. 
#and Algorithm 3.3 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear. 
#Requires (4m+2)n flops for m executions of the repeat-until loop if
#  sigma,Ms not given, otherwise 4mn flops.
#
#Reference implementation! Slow due to interpretation!
#

# written  03/03/07     S.M. Rump
#

  u = eps(S) / 2
  rm = realmin(S)
  if Ms == 0	                        # standard call
    n = length(p)                       # initialization
    mu = maximum(abs(p))                # abs(p_i) <= mu
    if n==0 || mu==0                    # no or only zero summands
      return rho, 0, p, 0, 0
    end
    Ms = NextPowerTwo(n+2)              # n+2 <= 2^M = Ms
    sigma = Ms * NextPowerTwo(mu)       # first extraction unit
    if ~isfinite(sigma)
      error("overflow occurred in Transform")
    end
  else                                  # sigma and Ms already known
    sigma = Ms * u * sigma
  end    
  phi = Ms * u                          # factor to decrease sigma    
  if kPhi == 0 Phi = Ms*Ms*u            # AccSum, TransformK, AccSumK, NearSum
  elseif kPhi == 1 Phi = Ms*u           # AccSign
  elseif kPhi == 2 Phi = 8*Ms*u         # AccSumHugeN
  end
  t = rho
  while true
    q = ( sigma + p ) - sigma           # [tau,p] = ExtractVector(sigma,p);
    tau = sum(q)                        # sum of leading terms
    p = p - q                           # remaining terms
    tau1 = t + tau                      # new approximation
    if ( abs(tau1) >= Phi * sigma ) || ( sigma <= rm )      
      tau2 = tau - ( tau1 - t )         # [tau1,tau2] = FastTwoSum(t,tau)
      return tau1, tau2, p, sigma, Ms
    end
    t = tau1                            # sum t+tau exact
    if t==0                             # accelerate case sum(p)=0
      tau1,tau2,p,sigma,Ms = Transform(p, 0, kPhi)  # sum of remainder part
      return tau1, tau2, p, sigma, Ms
    end
    sigma = phi * sigma                 # new extraction unit
  end
end
