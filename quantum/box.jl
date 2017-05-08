module Eigen

"""
Calculate the Energy discrete niveaus of a simply shaped Potential V.

V(x) = EL for x <= 0, V(x) = ER for x >= X, V(x) = EM for x in [0,X].

ER >= EL > EM.

Find eigenvalues E of Schrödinger Equation -hb²/(2m)*d²psi(x)/dx² + V(x) = E psi(x).

Discrete Eigenvalues E are expected with 0 < E <= min(EL, ER).
hb is Planck's constant divided by (2pi) (h-bar).

Restrictions: eM <= E <= min(EL, ER)
m mass of particle  ~9e-31 kg electron
X length parameter  ~5nm typical
E energy level      ~1e-18 J 

all SI units
hb = 6.626..e-34 J s

Zeros of determinant are the discrete eigenvalues
"""
function determinant(E, EL, EM, ER, m, X)
  hb = 6.62607004081e-34/(2pi) # J s
  x = sqrt((E-EM)*2m) * X / hb
  
  xL = sqrt((EL-EM)*2m) * X / hb
  xR = sqrt((ER-EM)*2m) * X / hb

  det_normal(x, xL, xR)
end

# normalized determinant. 0 <= x <= min(xL, xR)
function det_normal(x, xL, xR)
  xxL = sqrt(xL^2 - x^2)
  xxR = sqrt(xR^2 - x^2)
  sinc(x/pi) * (x^2 + sqrt(xxL * xxR)) + cos(x) * (xxR - xxL)
end

end # module
