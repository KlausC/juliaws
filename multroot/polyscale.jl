import Polynomials
using Polynomials

#
# given the coefficients of a polynomial p(x) = p(n)*x^n + + p(1)*x + p(0)
#
# determne a scaling factor for x and a linear factor t that the transformed
# polynomial q(z) = t * p(s*z) has coefficient p(n) = 1 and all |p(j)| <= 1
#
# Invariant:
# if q, shift, smax, u = polyscale(p) then q = polytransform(p, shift, smax) * u
# up to numerical precision

FloatOrComplex = Union{AbstractFloat, Complex, Complex{Float32}, Complex{BigFloat}}

function polyscale{T<:FloatOrComplex}(P::Poly{T})
	n = length(P.a) -1
    shift = - P.a[n] / n	
    P = polytransform(P, shift)	
	p = P.a
    p[n] = 0
    pni = 1 / p[n+1]
    smaxlog = maximum(log(abs(p[1:n] * pni)) ./ (n+1 - (1:n)))
    smax = exp(smaxlog)
    u = pni / exp(smaxlog * n)
    P = polytransform(Poly(p), 0.0, smax) * u
    P.a[n+1] = T(1)
    P, shift, smax, u
end
