function NextPowerTwo(p)
#NEXTPOWERTWO Smallest power of 2 not less than abs(p) for p~=0
#
#   L = NextPowerTwo(p)
#
#On return, abs(p) <= 2^L and L is smallest possible.
#
#Implements Algorithm 3.5 from
#  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
#    Faithful Rounding, to appear. 
#Requires 4 flops w/o scaling.
#
#Reference implementation! Slow due to interpretation!
#
#Improvement: in the following simpler version, the "if"-statement is omitted
#

# written  03/03/07     S.M. Rump
# modified 09/19/07     S.M. Rump  improved implementation: omit branch
#

  s, ex = frexp(float(p))
  ex + map( si -> si == 0.5 ? - 1 : 0, s)
end
