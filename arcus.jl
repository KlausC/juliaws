"""

  `module Arcus`

Module Arcus provides a data type Arc, representing points on the unit circle,
or an agular value in the range `(-π,π]`.
For a point given as `(sin(α), cos(α))`, the absolutely lower value of both is stored.
Quadrant information is stored by multiplying a corresponding power of two.
For the first quadrant with `α` in  `[-π/4,π/4]`, the Float64 value of `sin(α)` is
stored unchanged.

Usage:

  `using Arcus`

  `a = Arc(radians)`

  `b = Arc(s, c)`

Trigonometric functions and arithmetic operations are available.
Multiplication and division only with `Real` values.
"""
module Arcus

export Arc, Arcd, deg, rad, sincos
import Base: convert, show, zero, bits
import Base: sin, csc, cos, sec, tan, cot
import Base: (+), (-), (*), (/)

bitstype 64 Arc

const expomax = exponent(realmax(0.0)) ÷ 3 # use all available exponent space
const fmin = 2.0^-expomax
const f0 = 1.0                # quadrant 0 - stored values abs < f0
const f1 = 2.0^expomax        # quadrant 1 - f0 ≤ values < f1
const f2 = f1 * f1            # quadrant 2 - f1 ≤ values < f2
const f3 = f2 * f1            # quadrant 3 - f2 ≤ values < f3
const pi180 = pi / 180
const twopi = 2pi
const pi2 = pi / 2

"""

  `Arc(s::Real, c::Real) -> ::Arc`

Construct Arc from cartesian coordinates. `(s, c)` is forced to unit size.
"""
Arc(s::Real, c::Real) = reinterpret(Arc, sc_to_arc(s, c))

"""

  `Arc(rad::Real) -> ::Arc`

Construct Arc from radian value. If rad ∉ (-2π,2π] mind rounding errors.
"""
function Arc(r::Real)
  r = mod(r, twopi)
  if mod(r, pi2) == 0
    case = Int(fld(r, pi2))
    if case == 0
      Arc(0.0, 1.0)
    elseif case == 1
      Arc(1.0, 0.0)
    elseif case == 2
      Arc(0.0, -1.0)
    else
      Arc(-1.0, 0.0)
    end
  else
    s, c = sin(r), cos(r)
    Arc(s, c)
  end
end

"""

  `Arcd(degree::Real) -> ::Arc`

Construct Arc from degree value.
"""

# special exactness if d given as integer
function Arcd(d::Real)
  d = mod(d, 360)
  if mod(d, 90) == 0
    case = Int(fld(d, 90))
    if case == 0
      Arc(0.0, 1.0)
    elseif case == 1
      Arc(1.0, 0.0)
    elseif case == 2
      Arc(0.0, -1.0)
    else
      Arc(-1.0, 0.0)
    end
  else
    Arc(d * pi180)
  end
end

# working function for constructor
function sc_to_arc(s::Real, c::Real)
  iq = quadr(s, c)
  h = hypot(s, c)
  s /= h
  c /= h
  if iq == 0
    s
  elseif iq == 1
    ( abs(c) <= fmin ) && (c = fmin)
    c * f1
  elseif iq == 2
    ( abs(s) <= fmin ) && (s = fmin)
    s * f2
  else
    ( abs(c) <= fmin ) && (c = fmin)
    c * f3
  end
end

# Find quadrant number from s and c
function quadr(s::Real, c::Real)
  if abs(s) <= abs(c)
    ifelse( c > 0, 0, 2)
  else
    ifelse( s > 0, 1, 3)
  end
end

@inline sq(x::Real) = sqrt(one(x) - x * x)

"""

  `sincos(a::Arc) -> (s, c)::Tuple{Float64,Float64}`

Calculate both sine and cosine of `a`.
"""
function sincos(a::Arc)
  sc = reinterpret(Float64, a)
  sca = abs(sc)
  if sca < f0
    sc, sq(sc)
  elseif sca < f1
    sc /= f1
    ( abs(sc) <= fmin ) && return 1.0, 0.0
    sq(sc), sc
  elseif sca < f2
    ( abs(sc) <= fmin ) && return 0.0, -1.0
    sc /= f2
    sc, -sq(sc)
  else
    ( sca <= fmin ) && return -1.0, 0.0
    sc /= f3
    -sq(sc), sc
  end
end

"""

  `sin(a::Arc) -> s::Float64`

Calculate the sine of `a`.
"""
function sin(a::Arc)
  sc = reinterpret(Float64, a)
  sca = abs(sc)
  if sca <= f0
    sc
  elseif sca <= f1
    sc /= f1
    sq(sc)
  elseif sca <= f2
    sc /= f2
    ( abs(sc) <= fmin ) && return 0.0
    sc
  else
    sc /= f3
    -sq(sc)
  end
end

"""

  `cos(a::Arc) -> c::Float64`

Calculate the cosine of `a`.
"""
function cos(a::Arc)
  sc = reinterpret(Float64, a)
  sca = abs(sc)
  if sca <= f0
    sq(sc)
  elseif sca <= f1
    sc /= f1
    ( abs(sc) <= fmin ) && return 0.0
    sc
  elseif sca <= f2
    sc /= f2
    -sq(sc)
  else
    sc /= f3
    ( abs(sc) <= fmin ) && return 0.0
    sc
  end
end

"""

  `tan(a::Arc) -> c::Float64`

Calculate the tangent of `a`.
"""
function tan(a::Arc)
  s, c = sincos(a)
  s / c
end

"""

  `sec(a::Arc) -> c::Float64`

Calculate the secant of `a`.
"""
sec(a::Arc) = 1.0 / cos(a)
"""

  `csc(a::Arc) -> c::Float64`

Calculate the cosecant of `a`.
"""
csc(a::Arc) = 1.0 / sin(a)
"""

  `cot(a::Arc) -> c::Float64`

Calculate the cotangent of `a`.
"""
cot(a::Arc) = 1.0 / tan(a)

function convert(::Type{Float64}, a::Arc)
  s, c = sincos(a)
  atan2(s, c)
end

convert(::Type{Arc}, x::Real) = Arc(x)

"""

  `rad(a::Arc) -> r::Float`

Convert `a` to rad in `(-π,π]`.
"""
rad(a::Arc) = convert(Float64, a)
"""

  `deg(a::Arc) -> r::Float`

Convert `a` to degrees in `(180,180]`.
"""
deg(a::Arc) = rad(a) * 180 / pi

# Singular zero constant
const zero_arc = reinterpret(Arc, 0.0)
zero(::Type{Arc}) = zero_arc
zero(::Arc) = zero_arc

show(io::IO, a::Arc) =  print(io, "$(rad(a))r")

bits(a::Arc) = bits(reinterpret(Float64, a))

# Supported arithmetic functions
# plus
function +(a::Arc, b::Arc)
  sa, ca = sincos(a)
  sb, cb = sincos(b)
  Arc(sa*cb + sb*ca, ca*cb - sa*sb)
end

# minus
function -(a::Arc, b::Arc)
  sa, ca = sincos(a)
  sb, cb = sincos(b)
  Arc(sa*cb - sb*ca, ca*cb + sa*sb)
end

# unary -
function -(a::Arc)
  sa, ca = sincos(a)
  Arc(-sa, ca)
end

# times
function *(f::Integer, a::Arc)
  if abs(f) > 2
    Float64(f) * a
  elseif f == 2
    a + a
  elseif f == -2
    -( a + a )
  elseif f == -1
    -a
  elseif f == 1
    a
  elseif f == 0
    zero(Arc)
  end
end

*(f::Real, a::Arc) = Arc(rad(a) * f)
*(a::Arc, f::Real) = f * a

# diveded
/(a::Arc, f::Real) = (1.0 / f) * a

# Not required if not Arc <: AbstractFloat
# *(a::Arc, b::Arc) = error("invalid operation * with Arc")

end # module Arcus
