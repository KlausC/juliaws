
# Try to approximate erfc with the given function

f(t) = begin xt = x(t); log(erfcx(xt) / t) end
x(t) = 2 * (1/t - 1)
t(x) = 1 / (abs(x)/2 + 1)

include("approx.jl")

rat5 = approxfunction('a', degnom = 5)
erfcx_polynom(t, p) = t * exp(polynom(t, p))
erfc_rat5(t, p) = t * exp(-x(t)^2 + rat5(t, p))
erfcx_rat5(t, p) = t * exp(rat5(t, p))

ts = linspace(0, 1, 101)

erfcxt(t) = erfcx(x(t))
erfct(t) = erfc(x(t))


st, gp, fun, xopt, pp = approx(erfct, erfc_rat5, ts, [0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


