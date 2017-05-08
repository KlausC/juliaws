
function balance0{T<:Number}(A::AbstractArray{T,2}, tol::Float64 = sqrt(eps()))
  balance0(A, map( x -> x*x, A), tol)
end

function balance0{T<:Number}(A::AbstractArray{T,2}, A2::AbstractArray{T,2}, tol::Float64 = sqrt(eps()))
  m, n = size(A)
  println("Entering balance0 tol = ", tol)
  println("A = ", full(A))
  
 if m == 1
    xq = 1 ./ abs(vec(Array(full(sub(A,1,:)))))
    yq = ones(T,1)
 elseif n == 1
    xq = ones(T,1)
    yq = 1 ./ abs(vec(Array(full(sub(A,:,1)))))
 else
  xq = 1 ./ vec(maximum(abs(A), 1))
  yq = 1 ./ vec(maximum(abs(A), 2))
  x,y = xq, yq
  y = copy(yq)
  xv, yv = copy(xq), copy(yq)

  function dev(A, A2, x, y)
     na = length(A)
     a1 = sum(Diagonal(y) * A * Diagonal(x))
     a2 = sum(Diagonal(y .^ 2) * A2 * Diagonal(x .^ 2))
	 a2 - 2 * a1 + na
  end

  println("0: ", dev(A, A2, x, y))
  for l = 1:100

    for k = 1:n
      x1 = vecdot(sub(A,:,k), yq)
      x2 = vecdot(sub(A2,:,k), yq .^ 2)
      if x2 != 0 x[k] =  x1 / x2 end
    end

    for i = 1:m
      y1 = vecdot(sub(A,i,:), xq)
      y2 = vecdot(sub(A2,i,:), xq .^ 2)
      if y2 != 0 y[i] =  y1 / y2 end
    end
   
    xq = x 
    yq = y
    xb = sqrt(norm(xq) / norm(yq))
    xq /= xb
    yq *= xb
    println(l, ": ", dev(A, A2, xq, yq))
    if norm(qx-xv) < norm(xq) * tol && norm(yq-yv) <= norm(yq) * tol
      break
    end
    xv = copy(xq)
    yv = copy(yq)
  end
 end
 println("Leaving balance0")
 #println("x = ", xq)
 #println("y = ", yq)
 #println("A = ", Diagonal(yq) * full(A) * Diagonal(xq))
 xq, yq
end

function balance1{T<:Number}(A::AbstractArray{T,2}; tol::Float64 = sqrt(eps()), maxiter::Int = 100)
  m, n = size(A)
  maxx = 3
  println("Entering balance1 tol = ", tol)
  println("A = ", full(A))
  
  x = vec(one(T) ./ sum(A, 1))
  y = one(T) ./ ( A * x)

  function dev(A, x, y)
     sx = y' * A .* x' - size(A,1)
     sy = (A * x) .* y - size(A,2)
     norm(sx) + norm(sy)
  end

  println("0: ", dev(A, x, y))
  for l = 1:maxiter

    xinv = vec(y' * A)
    dxrel = norm(xinv .* x - m, Inf) / m
    if dxrel <= tol
      break
    end
    x = one(T) ./ xinv
    y = one(T) ./ (A * x) * n
   
    xb = sqrt(norm(x) / norm(y))
    x /= xb
    y *= xb
    #println(l, ": ", dev(A, x, y), " dx/x: ", dxrel)
  end
  #println("Leaving balance1")
  #println("x = ", xq)
  #println("y = ", yq)
  #println("A = ", Diagonal(yq) * full(A) * Diagonal(xq))
  x, y
end

function aitken(xacc)
  if size(xacc,2) < 3
    xacc[:,end]
  else
    xacc[:,end] - (xacc[:,end] - xacc[:,end-1]) .^ 2 ./ (xacc[:,end] - xacc[:,end-1] * 2 + xacc[:,end-2])
  end
end

  # construct summarizing smaller Matrix
  function summarize(A1, A2)
    m, n = size(A1)
    if m >= n
      m2 = (m+1)÷2
      B1 = spzeros(m2, n)
      B2 = spzeros(m2, n)
      for j = 1:m÷2
        B1[j,:] = abs(A1[j*2-1,:]) + abs(A1[j*2,:])
        B2[j,:] = abs(A2[j*2-1,:]) + abs(A2[j*2,:])
	  end
      if m % 2 == 1
        B1[m2,:] = abs(A1[m,:])
        B2[m2,:] = abs(A2[m,:])
      end
	  B1, B2
    else
      n2 = (n+1)÷2
      B1 = spzeros(m, n2)
      B2 = spzeros(m, n2)
      for k = 1:n÷2
        B1[:,k] = abs(A1[:,k*2-1]) + abs(A1[:,k*2])
        B2[:,k] = abs(A2[:,k*2-1]) + abs(A2[:,k*2])
	  end
      if n % 2 == 1
        B1[:,n2] = abs(A1[:,n])
        B2[:,n2] = abs(A2[:,n])
      end
	  B1, B2
    end
  end
 
  function balance_internal(A1, A2, tol)
    m, n = size(A1)
    xc, yc = 0.0, 0.0
	println("Entering balance_int m = ", m, " n = ", n)
    if min(m,n) <= 1
      xc, yc = balance1(A1)
    else
      B1, B2 = summarize(A1, A2)
      xb, yb = balance_internal(B1, B2, tol)
      xa = expa(xb, n)
      ya = expa(yb, m)
      xc, yc = balance1(Diagonal(ya) * A1 * Diagonal(xa)) 
    end
	println("Leaving balance_int")
    xc, yc
  end

  function expa(x, n)
    if length(x) == n
      x
    elseif n%2 == 0
      vec([x'; x'])
	else
	  [vec([x[1:end-1]'; x[1:end-1]']); x[end]]
    end
  end


function balance{T<:Number}(Ain::AbstractArray{T,2}, tol::Float64 = sqrt(eps()))

  colsum(A) = vec(sum(A .* A, 1) ./ sum(abs(A), 1))
  rowsum(A) = vec(sum(A .* A, 2) ./ sum(abs(A), 2))
  pc = sortperm(colsum(Ain))
  pr = sortperm(rowsum(Ain))
  A = sub(Ain, pr, pc) # A is sorted view on Ain - smallest first
  A2 = MapView(x->x.*x, A)

  balance_internal(A, A2, tol)
end

if ! isdefined(:MapView)
type MapView{T} <: AbstractArray{T}
  f::Function
  back::AbstractArray{T}
  MapView(f, A::AbstractArray{T}) = new(f, A)
end
MapView{T}(f::Function, A::AbstractArray{T}) = MapView{T}(f, A)
end
import Base: size, getindex, ndims, print_matrix, full, map
getindex(v::MapView, index...) = map(v.f, v.back[index...])
size(v::MapView) = size(v.back)
ndims(v::MapView) = ndims(v.back)
print_matrix(t1, v::MapView, t2...) = print_matrix(t1, map(v.f, v.back), t2...)
full(v::MapView) = full(map(v.f, v.back))
function map(g::Function, v::MapView) gf(x) = g(v.f(x))
  MapView(gf, v.back)
end

if ! isdefined(:VectorAccu)
type VectorAccu{T<:Number}
  a::AbstractArray{T,2}
  l::Int
  VectorAccu(m, len) = new(Array{Float64}(m, 0), len)
end
VectorAccu(m, n) = VectorAccu{Float64}(m, n)
end

function append!{T}(v::VectorAccu{T}, x::AbstractArray{T})
  if size(v.a,2) < v.l
    v.a = [v.a x]
  else
    v.a = [v.a[:,2:end] x]
  end
end

