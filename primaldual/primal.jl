
using MPS

function objective(p::MPS.LPData, x::Vector)
  size(x, 1) == p.n || error("Dimension mismatch")
  grad() = p.c' .* x'
  f(y::Vector) = grad() * y
  gradh() = p.A * Diagonal(x)
  h(y::Vector) = gradh() * y - p.b
  lambdah(y::Vector, lam::Vector) = abs(lam)' * abs(h(y))
  obj(y::Vector, lam::Vector) = lambdah(y, lam) + f(y)
  f, grad, h, gradh, lambdah, obj
end

using Base.SparseArrays.SPQR: Factorization, qmult, solve, QX, QTX, RETX_EQUALS_B, RTX_EQUALS_ETB
using Base.SparseArrays.CHOLMOD: Dense
function directions(QR::Factorization)
  function projgrad(g::Vector)
    Qtg = qmult(QTX, QR, Dense(g))
    prg = qmult(QX, QR, Qtg)
    g - prg
  end

  function restoration(h::Vector)
    rmth = solve(RTX_EQUALS_ETB, QR, Dense(h))
    d = qmult(QX, QR, rmth)
    convert(typeof(h), d)
  end
  
  projgrad, restoration
end
