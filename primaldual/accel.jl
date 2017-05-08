
if ! isdefined(:VectorAccu2)
type VectorAccu2{T<:Number}
  v::AbstractArray{T,2}
  w::AbstractArray{T,2}
  l::Int
  VectorAccu2(m, len) = new(Array{T}(m, 0), Array{T}(m,0), len)
end
VectorAccu2(T::Type, m, n) = VectorAccu2{T}(m, n)
end

function append!{T}(a::VectorAccu2{T}, x::AbstractArray{T}, y::AbstractArray{T})
  if size(a.v, 2) < a.l
    a.v = [a.v x]
    a.w = [a.w y]
  else
    a.v = [a.v[:,2:end] x]
    a.w = [a.w[:,2:end] y]
  end
  a
end

function pop!{T}(a::VectorAccu2{T})
  vn = a.v[:,end]
  wn = a.w[:,end]
  a.v = a.v[:,1:end-1]  
  a.w = a.w[:,1:end-1]  
  vn, wn
end

function differences{T}(a::VectorAccu2{T})
  if size(a.v, 2) <= 0
    error("accumulator is empty - no differences possible")
  end
  V = a.v[:,2:end] - a.v[:,1:end-1]
  W = a.w[:,2:end] - a.w[:,1:end-1]
  V, W
end

function accel{T}(a::VectorAccu2{T})
  Y, Z = differences(a)
  vn = a.v[:,end]
  wn = a.w[:,end]
  wv = wn - vn
  if norm(wv) > 0.0 && size(Y,2) > 0 && norm(Y[:,end]) > 0
    F = qrfact(Y)
    X = qrfact(eye(size(F,2)) - F \ Z)
    wn + Z * (X \ ( F \ wv))
  else
    wn
  end
end

#  Accu3 stores only differences. That will allow the superposition of values
#
if ! isdefined(:VectorAccu3)
type VectorAccu3{T<:Number}
  v::AbstractArray{T,2}
  w::AbstractArray{T,2}
  l::Int
  vn::AbstractVector{T}
  wn::AbstractVector{T}
  F::Factorization{T}
  function VectorAccu3(m, len)
    x = Vector{T}(m)
    new(Array{T}(m, 0), Array{T}(m,0), len, x, x)
  end
end
VectorAccu3(T::Type, m, n) = VectorAccu3{T}(m, n)
end

function append!{T}(a::VectorAccu3{T}, x::AbstractArray{T}, y::AbstractArray{T})
  n = size(a.v, 2)
  if  a.vn === a.wn && n == 0 # is it the first call to append! ?
    a.vn, a.wn = x, y
    a.F = qrfact(a.v) 
  else
    vn = x - a.vn
    wn = y - a.wn
    rn = (vn' * a.F[:Q])[1:size(a.F,2)]'
    nrn = vecnorm(rn)
    nvn = vecnorm(vn)
    # check if new vector is too close to span of previous vectors
    if n > 0 && (n >= a.l || nrn >= 0.8 * nvn || nvn == 0.0 )
      # spread new entry over all existing ones
      a.v = a.v + vn * rn
      a.w = a.w + wn * rn
    else
      a.v = [a.v vn]
      a.w = [a.w wn]
    end
    a.vn = x
    a.wn = y
    a.F = qrfact(a.v)
  end
  a
end

function accel{T}(a::VectorAccu3{T})
  Y, Z = a.v, a.w
  wv = a.wn - a.vn
  if norm(wv) > 0.0 && size(Y,2) > 0 && norm(Y[:,end]) > 0
    F = qrfact(Y)
    X = qrfact(eye(size(F,2)) - F \ Z)
    a.wn + Z * (X \ ( F \ wv))
  else
    a.wn
  end
end
