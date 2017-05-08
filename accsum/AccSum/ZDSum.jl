function ZDSum{S}(v::AbstractArray{S,})
#ZDSUM        Computes approximation of sum(p) using Zielke/Drygalla's approach
#
#   y = ZDsum(p)
#
#On return, res is an approximation of sum(p) of the column vector p
#  provided no over- or underflow occurred. 
#
#Taken from Algorithm accdot on page 30 in
#  G. Zielke, V. Drygalla: Genaue Lösung linearer Gleichungssysteme, 
#    GAMM Mitteilungen 26(1-2), p. 8-107, 2003.
#Requires (6m+3)n+7k flops for m executions of while-loop and
#  length(s)=k for intermediate vector s.
#
# written  03/03/07     S.M. Rump
#

  n = length(v)
  if S == Float64
    C1 = 54;
    C2 = -1023;
  else
    C1 = 25;
    C2 = -127;
  end
  ma = maximum(abs(v))
  emax = exponent(ma) # [mmax, emax] = log2(ma)
  q = 2^emax
  k = floor(C1 - log2(n))
  p = 2^k
  i = 0
  s = zeros(S,0)
  while any(vi -> vi != 0, v) && ( q / p > 2 ^ C2 )
    i = i+1
    q = q / p
    g = round(Int, v / q, RoundToZero)	# fix(v/q)
    push!(s, sum(g))
    v = v - g * q
  end
  i = i+1
  push!(s, sum(v) * p / q)
  ue = 0
  for j = i:-1:1
    t = s[j] + ue;
    ue = floor(t / p)
    s[j] = t - ue * p
  end
  y = ue
  for j = 1:i
    y = s[j] + y * p
  end
  y = y * q / p
  y
end
