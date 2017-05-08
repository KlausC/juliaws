function FastPrecSum{S<:Number}(p::AbstractArray{S,}, K::Integer)
#FastPrecSum  Ultimately fast and accurate summation as if computed in K-fold precision
#
#   res = fastprecsum(p)
#
#For real or complex input vector, dense or sparse, the result res is
#sum(p_i) as if computed in K-fold precision. Input vector p must not be of type intval.
#
#Maximum number of nonzero elements per sum is limited to 67,108,862 in 
#double precision, which seems sufficient for Matlab.
#
#Implements new algorithm in
#  S.M. Rump: Ultimately Fast Accurate Summation, submitted for publication, 2008.
#
#CAUTION: !!! THIS IMPLEMENTATION SUFFERS SEVERELY FROM INTERPRETATION OVERHEAD !!!
#!!! IT IS INCLUDED TO SHOW THE PRINCIPLES OF THE NEW METHOD !!!
#!!! DO NOT USE FOR LARGE DIMENSIONS !!!
#

# written  10/18/08     S.M. Rump
#

  res = 0
  m = 0
  if isempty(p)
    return res
  end

  # check size
  if length(size(p)) > 2
    error("fastprecsum not defined for multi-dimensional arrays.")
  end
  
  # check improper input
  if any(isnan(p)) || any(isinf(p))
    res = NaN
    return res
  end

  # take care of complex input
  if !isreal(p)
    resreal,exactreal = fastprecsum(real(p))
    resimag,exactimag = fastprecsum(imag(p))
    exact = exactreal && exactimag
    res = resreal + resimag*im
    return res
  end

  # input real, compute sum(p)
  e = 1e-30
  if 1+e == 1-e                           # fast check for rounding to nearest
    rndold = 0
  else
    rndold = getround ## TODO
    setround(0)
  end

  if issparse(p)
    n = nnz(p)                         # initialization
  else
    n = length(p)
  end
  
  # initialize constants depending on precision
  if S == Float32
    eps = 2^(-24)
    eta = 2^(-149)
  else
    eps = 2^(-53)
    eta = 2^(-1074)
  end
  
  # check dimension
  if ((2n+4) * n + 6) * eps > 1
    error("dimension too large for fastprecsum")
  end
  
  # initialize constants
  ExactFlag = false
  c1 = 1 - n * eps  
  c2 = 1 - (3n + 1) * eps
  c3 = 2eps
  c4 = 1 - eps  
  c5 = 2n * (n+2) * eps  
  c6 = 1 - 5eps  
  c7 = (1.5 + 4eps) * (n * eps)
  c8 = 2n * eps
  c9 = eta / eps
  sigma0 = zeros(1,K)
  Phi = sigma0
  
  T = sum(abs(p)) / c1                  # sum(abs(p)) <= T
  if T <= c9                            # no rounding error
    res = sum(p)
    if rndold setround(rndold) end
    return res, m
  end
  
  # intialize constants
  for m = 1:K-1
    sigma0[m] = 2T / c2                 # array T not needed
    q = sigma0[m] / c3
    u = abs(q / c4 - q)                 # u = ufp(sigma0)
    Phi[m] = ( c5 * u ) / c6
    T = min( c7*sigma0[m] , c8*u )      # sum(abs(p)) <= T
    if 4T <= c9
      K = m + 1
      ExactFlag = true
	end
  end
  sigma0p = sigma0
  
  # extraction loop
  e = 0;
  for i = 1:n
    ps = p[i]
    for m = 1:K-1
      sigmap = sigma0p[m] + ps
      q = sigmap - sigma0p[m]
      ps = ps - q
      sigma0p[m] = sigmap
    end
    e = e + ps
  end
  
  if ExactFlag                          # sum exact
    return e
  end
  
  # final result
  t = 0
  for m = 1:K-1
    tau = sigma0p[m] - sigma0[m]        # exact result
    tm = t + tau
    if abs(tm) >= Phi[m]                # surely faithful rounding
      tau2 = (t - tm) + tau             # [tm,tau2] = FastTwoSum(t,tau)
      if m == K-1
        res = tm + (tau2 + e)
      else                              # m <= K-2
        taum1 = sigma0p[m+1] - sigma0[m+1]
        tau3 = tau2 + taum1
        tau4 = (tau2 - tau3) + taum1    # [tau3,tau4] = FastTwoSum(tau2,taum1)
        if m == K-2
          res = tm + (tau3 + (tau4 + e))
        else
          res = tm + (tau3 + (tau4+ (sigma0p[m+2] - sigma0[m+2])))
        end
      end
      return res
    end
    t = tm
  end
  
  res = tm + e                         # maybe not faithful
  res
end
