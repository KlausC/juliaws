function PriestSum(p)
#PRIESTSUM    Priest's doubly compensated summation
#
#   res = PriestSum(p)
#
#Approximation of sum(p) almost accurate to two units of the last place.
#
#Follows D.M. Priest: On properties of floating point arithmetics: Numerical stability 
#           and the cost of accurate computations, Mathematics Department, University 
#           of California at Berkeley, CA, USA (1992).
#Requires n log(n) flops.
#
# written  03/03/07     S.M. Rump
#

  index = sortperm(abs(p), rev = true)
  p = p[index]
  s = p[1]
  c = 0
  for i = 2:length(p)
    y,u = FastTwoSum(c, p[i])
    t,v = FastTwoSum(s, y)
    z = u + v
    s,c = FastTwoSum(t, z)
  end
  s
end
