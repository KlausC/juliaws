module Polygon

using Base.Random

"""
Create a random-walk polygon, returning to initial point.
"""
function closed_walk{F<:AbstractFloat}(n::Int, T::Type{F} = Float64; rng = GLOBAL_RNG)
  x = randn(rng, T, n)
  y = randn(rng, T, n)
  xya = [hypot(p[1],p[2]) for p in zip(x,y)]
  x ./= xya
  y ./= xya
  x -= mean(x)
  y -= mean(y)
  x = cumsum(x)
  y = cumsum(y)
  x -= mean(x)
  y -= mean(y)
  x, y
end

"""
Create a random-walk polygon, returning to initial point, which does not operlap.
"""
function closed_walk_no{F<:AbstractFloat}(n::Int, T::Type{F} = Float64; utol = 1.5, ftol = 0.4, rng = GLOBAL_RNG)
  x = randn(rng, T, n)
  y = randn(rng, T, n)
  xya = [hypot(p[1],p[2]) for p in zip(x,y)]
  x ./= xya
  y ./= xya
  rx = Vector{T}(n)
  ry = Vector{T}(n)
  rx[1] = x[1]; ry[1] = y[1]
  for k = 2:n
    ax = rx[k-1]; ay = ry[k-1]
    fact = zero(T)
    while fact < ftol
      bx = ax + x[k]; by = ay + y[k]
      uminpos = Inf
      umaxneg = -Inf
      for j = 2:k-2
        cx = rx[j-1]; cy = ry[j-1]
        dx = rx[j]; dy = ry[j]
        u, v = intersect_point(ax,ay,bx,by,cx,cy,dx,dy)
        if -0.5 <= v <= 1.5
          println("k = $k j = $j u = $u")
          if u >= 0.0
            uminpos = min(uminpos, u)
          else
            umaxneg = max(umaxneg, u)
          end
        end
      end
      println("k = $k u-Interval [$umaxneg,$uminpos]")
      uabs = max(uminpos, -umaxneg)
      fact = one(T)
      if uabs < utol
        if uminpos >= -umaxneg
          fact = uminpos / utol
        else
          fact = umaxneg / utol
        end
      else
        if uminpos < utol
          fact = -fact
        end
      end
      if fact < ftol
        x[k], y[k] = randn(rng, T, 2)
        xy = hypot(x[k], y[k])
        x[k] /= xy
        y[k] /= xy
      end
    end
    println("k = $k fact = $fact")
    rx[k] = ax + x[k] * fact
    ry[k] = ay + y[k] * fact
  end
  rx, ry
end

"""
Create regular polygon with radius 1, random phase factor
"""
function regular{F<:AbstractFloat}(n::Int, T::Type{F} = Float64; rng = GLOBAL_RNG)
  ph = rand(rng) * 2pi
  r = 0:(n-1)
  x = map(k -> T(cos(2pi*k/n+ph)), r) 
  y = map(k -> T(sin(2pi*k/n+ph)), r) 
  x, y
end

"""
Given 2-dim points a,b,c,d find u and v so a + u*b = c + v*d
If no such u, v exist set u = v = infinite.
If infinitely many such u,v exist, select one with v = 0
"""
function intersect_point(ax::Real,ay::Real,bx::Real,by::Real,cx::Real,cy::Real,dx::Real,dy::Real)
  bx -= ax; by -= ay
  dx -= cx; dy -= cy
  cx -= ax; cy -= ay
  det = by * dx - bx * dy
  if det != 0
    u = (dx * cy - dy * cx) / det
    v = (bx * cy - by * cx) / det
  elseif bx * cy == by * cx
    v = zero(det)
    if abs(bx) >= abs(by)
      u = cx / bx
    else
      u = cy / by
    end
  else
    u, v, _ = promote(Inf, Inf, det)
  end
  u, v
end

end #module
