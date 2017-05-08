function FastAccSum{S<:Number}(p::AbstractArray{S,})
#FastAccSum   Ultimately fast and accurate summation with faithful rounding
#
#   res, m = fastaccsum(p)
#
#For real or complex input vector, dense or sparse, the result res is
#sum(p_i) faithfully rounded. Input vector p must not be of type intval.
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

# written  08/28/08     S.M. Rump
# modified 09/28/08     S.M. Rump  check for rounding to nearest improved
#

  res = 0
  m = 0
  if isempty(p)
    return
  end

  # check size
  if length(size(p)) > 2
    error("fastaccsum not defined for multi-dimensional arrays.")
  end
  
  # check improper input
  if any(isnan(p)) | any(isinf(p))
    res = NaN;
    return
  end

  # take care of complex input
  if !isreal(p)
    resreal,exactreal = fastaccsum(real(p))
    resimag,exactimag = fastaccsum(imag(p))
    exact = exactreal && exactimag
    res = resreal + resimag*im
    return res, exact
  end

  # input real, compute sum(p)
  e = 1e-30
  if 1+e == 1-e                           # fast check for rounding to nearest
    rndold = 0
  else
    rndold = getround
    setround(0) ## TODO
  end

  if issparse(p)
    n = nnz(p)                            # initialization
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
  if ((2*n+4) * n + 6) * eps > 1
    error("dimension too large for fastaccsum")
  end
  
  # initialize constants
  c1 = 1 - n * eps
  c2 = 1 - (3n+1) * eps
  c3 = 2eps
  c4 = 1 - eps
  c5 = 2n * (n+2) * eps
  c6 = 1 - 5eps
  c7 = (1.5 + 4eps) * (n * eps)
  c8 = 2n * eps
  c9 = eta / eps
  m = 0
  
  T = sum(abs(p)) / c1                  # sum(abs(p)) <= T
  if T <= c8                            # no rounding error
    res = sum(p)
    if rndold  setround(rndold)  end
    return res, m
  end
  tp = 0
  while true
    m = m + 1
    sigma0 = (2T) / c2
    P = cumsum([sigma0 p])              # [sigma_n,p] = ExtractVectorNew(sigma0,p)     
    q = P[2:n+1] - P[1:n]
    p = p - q                           # extracted vector
    tau = P[n+1] - sigma0               # tau = sigma_n-sigma0 exact
    t = tp
    tp = t + tau                        # s = t + tau + sum(p)
    if tp == 0                          # check for zero t+tau
      res,M = FastAccSum(p[find(p)])    # recursive call, zeros eliminated
      m = m + M
      if rndold setround(rndold) end
      return res, m
    end
    q = sigma0 / c3
    u = abs(q/c4 - q)                   # u = ufp(sigma0)
    Phi = ( c5 * u ) / c6
    T = min( c7 * sigma0 , c8 * u )     # sum(abs(p)) <= T
    if ( abs(tp) >= Phi ) || ( 4T <= c9 )
      tau2 = (t-tp) + tau               # [tp,tau2] = FastTwoSum(t,tau)
      res = tp + ( tau2 + sum(p) )      # faithful.y rounded result
      if rndold setround(rndold) end
      return res, m
    end
  end
end
