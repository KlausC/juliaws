function PrecSum{S}(p::AbstractArray{S,}, K::Integer = 2)
#PRECSUM      Fast summation better than K-fold precision
#
#   res = PrecSum(p,K)
#
#On return, either res is a faithful rounding of s:=sum(p), or
#
#   abs(res-s) < eps^K sum(abs(p)) .
#
#The default is K=2.
#
#The statements are also true in the presence of underflow. 
#The input vector p may be single or double precision.
#
#Implements Algorithm 4.1 from
#  S.M. Rump, T. Ogita, S. Oishi: Fast High Precision Summation, to appear. 
#Requires (4L+3)n flops for L as computed in the algorithm. For double precsion
#and the default K=2, we have L=2 or L=3. Note that an implementation in C should
#use the modification as in Algorithm 5.2 in the cited paper to avoid extra memory.
#
# written  09/29/07     S.M. Rump
#

  nmax = 2^div(Base.significand_bits(S)+1,2) - 2
  rm = realmin(S)
  n = length(p)
  if n > nmax
    error("maximum length of input vector $(nmax)")
  end
  
  u = 0.5*eps(S)
  mu = sum(abs(p))/(1-2*n*u)
  if mu==0 res=0; return res end
  Ms = NextPowerTwo(n+2)                # n+2 <= 2^M = Ms
  M = log2(Ms)
  L = ceil( ( K*log2(u) - 2 ) / ( log2(u) + M ) ) - 1
  sigma = zeros(L+1)                    # sigma_k in the paper is sigma(k+1)
  tau = zeros(L)                        # hosts sum of leading terms
  sigma[1] = NextPowerTwo(mu)           # first extraction unit
  if !isfinite(sigma)
    error("overflow occurred in PrecSum")
  end
  
  phi = Ms * u                          # factor to decrease sigma
  for k = 1:L                           # compute sigma(i)
    if sigma[k] > rm
      sigma[k+1] = phi * sigma[k]
    else
      L = k-1; break
    end
  end
  
  if L == 0                             # sigma_0 in underflow range
    res = sum(p)
    return res
  end
  
  for k = 1:L                           # main loop: vertical version
    q = ( sigma[k] + p ) - sigma[k]     # [tau,p] = ExtractVector(sigma_{k-1},p);
    tau[k] = sum(q)                     # sum of leading terms
    p = p - q                           # remaining terms
  end
  
  pi = tau[1]; e = 0;
  for k = 2:L                           # compute final result
    x = pi;
    pi = x + tau[k]                     # [pi,q]=FastTwoSum(pi,tau(k))
    q = tau[k] - ( pi - x )
    e = e + q                           # fl-pt sum of errors
  end
  res = pi + ( e + sum(p) )             # final result
  res
end
