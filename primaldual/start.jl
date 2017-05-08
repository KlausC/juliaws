
cd(ENV["HOME"] * "/juliaws/primaldual")

include("setup.jl")
include("projgrad.jl")
include("conjgrad.jl")

using Optim
using PyPlot

pr = problem("maros-r7")
m, n = size(pr.A)

f, h, Df, Dh, D2L = functions(pr)
g, Dg, D2g = penalty_functions()

function setx(x::Vector{Float64})

  B, T, Tl, d, lambda, TLT, TLTm1, TTLTT, fC = dirfuncs(x)

  dx = d(h(x))
  gx = TTLTT(Dg(x))
  nx = TTLTT(Df(x))

  duvw{U <: Real, V <: Real, W <: Real}(u::U, v::V, w::W) = dx * u + gx * v + nx * w
  xuvw(u, v, w) = x - duvw(u, v, w)

  fval(f, u, v, w) = f(xuvw(u, v, w))

  function stepmax{T <: AbstractVector{Float64}}(x, d::T)
    mx = 1.0 / maximum(d ./ x)
    mx >= 0 ? mx : Inf
  end
  stepmax(d) = stepmax(x, d)

  B, T, Tl, d, lambda, TLT, TLTm1, TTLTT, fC, dx, gx, nx, duvw, xuvw, fval, stepmax
end

x = ones(n)
B, T, Tl, d, lambda, TLT, TLTm1, TTLTT, fC, dx, gx, nx, duvw, xuvw, fval, stepmax = setx(x)

history = Dict{Int, Tuple}()

function iterate(x, alpha, stp)
  B, T, Tl, d, lambda, TLT, TLTm1, TTLTT, PRE, dx, gx, nx, duvw, xuvw, fval, stepmax = setx(x)
  tx = nx + alpha * gx
  res = optimize(t -> f(x-t*tx) + alpha * g(x-t*tx), 0.0, stepmax(x, tx))
  x = x - res.minimum * tx
  history[stp] = (extrema(x), f(x), g(x), res.f_minimum, alpha, res.minimum, x)
  x
end

