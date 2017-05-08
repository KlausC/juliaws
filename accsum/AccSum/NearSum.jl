function NearSum{S}(p::AbstractArray{S,})
#NEARSUM      Computes the rounding to nearest result of sum(p)
#
#   resN = NearSum(p)
#
#On return, resN is sum(p) rounded to nearest, also in the presence
#  of underflow. Input vector p may be single or double precision.
#
#Implements Algorithm 7.3 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
#    Sign, K-fold Faithful and Rounding to Nearest, to appear. 
#Requires (4m+4k+4)n flops for m,k executions of the repeat-until loops
#  in the calls of TransformK, respectively.
#
#Reference implementation! Slow due to interpretation!
#

# written  03/03/07     S.M. Rump
#

  nmax = 2^div(Base.significand_bits(S)+1,2) - 2
  if length(p) > nmax
    error("maximum length of input vector for NearSum $(nmax) .")
  end

  eta = nextfloat(0*p(1))
  tau1,tau2,p,sigma,Ms = Transform(p,0)
  tau2s = tau2 + sum(p)
  res,delta = FastTwoSum(tau1, tau2s)       # res+delta = tau1+tau2s
  if delta == 0                             # res = fl(s)
    return res
  end
  R = tau2 - ( res - tau1 )                 # s-res = R+sum(p)
  if delta < 0                              # fl(s) in {pred(res),res}
    gamma = prevfloat(res) - res            # res+gamma=pred(res)
    if gamma == -eta                        # s = res
      return res
    end
    deltas = gamma / 2                      # mu := res+deltas = M^-(res)
    deltass = TransformK(p,R-deltas,sigma,Ms)  # s-mu = R-deltas+sum(p)
    if deltass > 0                          # s > M^-(res)
      resN = res
    elseif deltass < 0                      # s < M^-(res)
      resN = prevfloat(res)
    else                                    # s = M^-(res)
      resN = res + deltas
    end
  else                                      # fl(s) in {res,succ(res)}
    gamma = nextfloat(res) - res            # res+gamma=succ(res)
    if gamma == eta                         # s = res
      return res
    end
    deltas = gamma / 2                      # mu := res+deltas = M^+(res)
    deltass = TransformK(p,R-deltas,sigma,Ms)  # s-mu = R-deltas+sum(p)
    if deltass > 0                          # s > M^+(res)
      resN = nextfloat(res)
    elseif deltass < 0                      # s < M^+(res)
      resN = res
    else                                    # s = M^+(res)
      resN = res + deltas
    end
  end
  resN
end
