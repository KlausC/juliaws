
include("fejerclenshaw2.jl")

function  makesums(n::Int)
  
  function sumprep(f, a::Real, b::Real)

    f2(t) = f(a*(1-t) + b*(1+t))
    fac = b - a
    f2, fac
  end

  function sum(f, a, b, xw)
    x = xw[:,1]
    w = xw[:,2]
    f2, fac = sumprep(f, a, b)
    fac * dot(w, map(f2, x))
  end

  wf1, wf2, wcc = fejer2(n)

  f1(f, a, b) = sum(f, a, b, wf1)
  f2(f, a, b) = sum(f, a, b, wf2)
  f3(f, a, b) = sum(f, a, b, wcc)
  f1, f2, f3
end


