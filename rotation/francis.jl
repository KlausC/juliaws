"""
Calculation of eigenvalues using Francis Iteration - QR
"""
module Francis

import Base: A_mul_Bc!, A_mul_B!, (*)


typealias AbstractM Union{AbstractMatrix, LinAlg.AbstractRotation}

"""
One step of Francis transformation Transform hessenberg matrix
H is in upper Hessenberg form.
shifts is a vector of shift values.
With each shift value, its conjugate complex value is implicitly used as shift
All shift values should be distinct.
"""
function transform_Hess!{T<:Real, S<:Union{Real,Complex}}(A::AbstractMatrix{T}, ilo::Integer, ihi::Integer, Q::AbstractM, s::Vector{S}, maxchase::Integer)

  !isa(Q, AbstractMatrix) || size(A, 2) == size(Q, 2) || error("A and Q have incompatible sizes")
  n, m = size(A)
  z = zero(T)
  n == m || error("require square matrix")
  1 <= ilo <= ihi <= n || error("ilo = $ilo, ihi = $ihi not in ($n, $m) Matrix)")

  m = counteigs(s)
  #ilo, ihi = deflate_subdiagonal!(A, ilo, ihi)
  #println("def_sub ilo = $ilo ihi = $ihi")
  
  # supress further calculations if A is already quasi-diagonal
  isize = ihi - ilo + 1
  if isize <= 1 || isize == 2 && discriminant(A, ilo, ihi) < z
    m = 0
    maxchase = 0
  end
  
  # Any shifts provided
  if m > 0

    println("before initial chase m = $m"); display(A)

    # First step: shift bulges blocking the undisturbed insertion
    if isize > 4
      i0, m0 = lastbulge(A, ilo, ihi, m + 1)
      if i0 > 0
        chase!(A, ilo, ihi, Q, m + 2 - i0)
      end
    end

    println("after initial chase m = $m"); display(A)
    # Second step: compute  prod(A-s(k))e1
    pe = zeros(T, n)
    pe[ilo] = one(T)
    k = j = 1
    while j <= m
      sr = real(s[k])
      si = imag(s[k])
      if si == 0
        pe = A * pe - sr * pe
        j += 1
      else
        pe = A * ( A * pe - 2sr * pe) + ( hypot(sr, si) ^ 2 ) * pe
        if k < length(s) && imag(s[k+1]) == -si
          k += 1
        end
        j += 2
      end
      k += 1
    end

    # Third step: set in new upper left bulge
    m2 = findlast(x -> x != 0, pe)
    for k = m2:-1:ilo+1
      G, r = givens(pe, k-1, k)
      pe[k-1] = r
      pe[k] = 0
      A_mul_B!(G, A)
      A_mul_Bc!(A, G)
      A_mul_Bc!(Q, G)
    end
  end
  # println("m = $m") 
  # display(A)

  # Forth step: chase all remaining bulks, oldest first
  if maxchase > 0
    display(A)
    println("chase rem")
    chase!(A, ilo, ihi, Q, maxchase)
  end
  ilo, ihi, Q
end

function transform_Hess!{T<:Real, S<:Union{Real,Complex}}(A::AbstractMatrix{T}, Q::AbstractM, s::Vector{S}, maxchase::Integer)

  transform_Hess!(A, 1, size(A, 1), Q, s, maxchase)
end

function transform_Hess!{T<:Real, S<:Union{Real,Complex}}(A::AbstractMatrix{T}, s::Vector{S}, maxchase::Integer)

  #Q = LinAlg.Rotation(LinAlg.Givens{T}[])
  Q = eye(T, size(A, 2))
  transform_Hess!(A, Q, s, maxchase)
end

"""
Count number of eigenvalues.
Each non-real eigenvalue is counted twice, if not paired with its conjugate.
"""
counteigs{T<:Real}(s::Vector{T}) = length(s)
function counteigs{T}(s::Vector{Complex{T}})
  z = zero(T)
  m = length(s)
  count = 0
  k = 1
  while k <= m
    sk = s[k]
    si = imag(sk)
    if si == z
      count += 1
    else
      if k < m && imag(s[k+1]) == -si
        k += 1
      end
      count += 2
    end
    k += 1
  end
  count
end

"""
Chase all bulges towards right lower corner
"""
function chase!(A, ilo::Integer, ihi::Integer, Q, maxchase)
  n1, n2 = size(A)
  n2 = min(ihi, n2)
  n1 = min(ihi, n1)
  m = i0 = n2 - 1 # column to start with

  while i0 >= ilo && m > 0
    i0, m = lastbulge(A, ilo, ihi, i0 - 1)
    println("i0 = $i0, m = $m, ilo = $ilo  ihi = $ihi")
    if m > 0
      for  i = i0:min(i0 + maxchase - 1, n2 - 2)
        for k = min(i+m+1,n1):-1:i+2
          if A[k,i] != 0
            G, r = givens(A, k-1, k, i)
            A_mul_B!(G, A)
            A_mul_Bc!(A, G)
            A_mul_Bc!(Q, G)
            A[k-1,i] = r
            A[k,i] = 0
          end
          println("after i = $i k = $k")
          display(A)
        end
      end
    end
  end
  A, Q
end

"""
Find last bulge start column and size
"""
function lastbulge{T<:Number}(A::AbstractMatrix{T}, ilo::Integer, ihi::Integer, i::Integer)
  n, m = size(A)
  n == m || error("need square matrix")
  #i <= n - 1 || error("column number must be less than n")

  maxpos = 2
  while i >= ilo
    m = findlast(A[i+2:ihi,i])
    if m == 0 
      if maxpos > 2
        break
      end
    else
      maxpos = max(maxpos, m + i + 1)
    end
    i -= 1
  end
  m = maxpos - i - 2
  m == 0 ? 0 : i + 1, m
end

"""
Find estimations for eigenvalues from lower right.
"""
function estimations!(A::AbstractMatrix, ilo::Integer, ihi::Integer, Q::AbstractM, isize::Integer)
  n, m = size(A)
  n == m || error("square Matrix required")
  n = min(n, ihi)
  1 <= isize <= n || error("submatrix not ok")
  ni = max(ilo-1, n - isize)
  isize = ihi - ilo
  if isize == 1
    return A
  elseif isize == 2
    if discriminant(A, ihi-1, ihi) < 0
      return A
    end
  end
  F = schurfact(view(A, ni+1:n, ni+1:n)) # does not modify A
  transform!(A, ilo, ihi, Q, F[:Z], ni, ni)
end

"""
Discriminant of 2 x 2 submatrix
"""
@inline discriminant(A, i, j) = ((A[i,i] - A[j,j]) / 2 ) ^2 + A[i,j] * A[j,i]


function transform!(A::AbstractMatrix, ilo::Integer, ihi::Integer, Q::AbstractM, FZ::AbstractMatrix, ni::Integer, mi::Integer)
  n, m = size(A)
  n = min(n, ihi)
  A[ni+1:n, mi:m] = FZ * A[ni+1:n, mi:m]
  A[ilo:n,mi+1:n] *= FZ' 
  Q[ilo:n,mi+1:n] *= FZ'
  A
end

"""
Deflate elements in first subdiagonal.
Return new lower and upper indices.
"""
function deflate_subdiagonal!(A::AbstractMatrix, ilo::Integer, ihi::Integer)
  z = zero(eltype(A))
  while ihi > ilo && delfation_criterion(A[ihi,ihi-1], A[ihi-1,ihi-1], A[ihi,ihi])
    A[ihi,ihi-1] = z
    ihi -= 1
  end
  while ilo < ihi && delfation_criterion(A[ilo+1,ilo], A[ilo+1,ilo+1], A[ilo,ilo])
    A[ilo+1,ilo] = z
    ilo += 1
  end
  ilo, ihi
end

"""
Early aggressive defaltion step without re-ordering
Return new upper index.
"""
function deflate!(A::AbstractMatrix, ilo::Integer, ihi::Integer, isize::Integer)
  z = zero(eltype(A))
  n, m = size(A)
  n, m = min(n, ihi), min(m, ihi)
  mi, ni = max(1, m - isize), max(1, n - isize)
  v = A[ni+1:n, mi]
  for k = n:-1:ni+1
    if deflation_criterion(A[k,mi], A[k,k], A[mi,mi])
      A[k,mi] = z
      ihi -= 1
    else
      break
    end
  end
  ihi
end

"""
deflation criterion.
"""
@inline deflation_criterion{T}(sub::T, da::T, db::T) = abs(sub) <= (abs(da) + abs(db) ) * eps(T)

"""
Repeated aggressive deflation steps re-ordering Eigenvalues
"""
function reorder!(A::AbstractMatrix, ilo::Integer, ihi::Integer, Q::AbstractM, isize::Integer)
  n = size(A, 1)
  n = min(n, ihi)
  ni = n - isize + 1
  to = ihi
  while ( (from, n1, to, n2) = lastpairs(A, ni, to); n1 ) > 0
    swap!(A, ilo, ihi, Q, from, n1, to, n2)
    ni = from + n2
    to = ni - 1
  end
  ihi = deflate!(A, ilo, ihi, isize)
  ni, ihi
end

"""
Find last pair of eigenvalue blocks in Schur matrix diagonal.
"""
function lastpairs(A::AbstractMatrix, ilo::Integer, ihi::Integer)
  z = zero(eltype(A))
  n2 = ihi > ilo && A[ihi,ihi-1] != z ? 2 : 1
  ihi -= n2
  n1 = ihi > ilo && A[ihi,ihi-1] != z ? 2 : ihi >= ilo ? 1 : 0
  to = ihi - n1 + 1
  to, n1, ihi + 1, n2
end

"""
Decrease index i, if subdiagonal element left of i,i is not zero
"""
@inline function ibottom(A::AbstractMatrix, ilo::Integer, iup::Integer)
  iup <= ilo || A[iup,iup-1] == 0 ? iup : iup - 1
end

"""
Calculate start index and sizes of diagonal blocks.
"""
function index_sizes(A::AbstractMatrix, ilo::Integer, iup::Integer)
  i2 = ibottom(A, ilo, iup)
  n2 = iup - i2 + 1
  i1 = ibottom(A, ilo, i2 - 1)
  n1 = i2 - i1
  return i1 < ilo ? (0, 0, 0) : (i1, n1, n2)
end

"""
Perform swap and check if spike vector is improved by the swap.
"""
function testswap(A::AbstractMatrix, spike::AbstractVector, ilo::Integer, iup::Integer)
  i1, n1, n2 = index_sizes(A, ilo, iup)
  r = i1:iup
  is_swap_ok(A[r,r], spike[r], n1, n2)
end
"""
A is a n1+n2:n1+n2 matrix, spike a n1+n2 vector.
"""
function is_swap_ok{T<:Real}(A::AbstractMatrix{T}, spike::AbstractVector{T}, n1::Integer, n2::Integer)

  ze, on = zero(T), one(T)
  g(i1::Integer, i2::Integer) = givens(ze, on, i1, i2)[1]
  n = n1 + n2 
  # source ranges
  spnorm1 = norm(spike[n1+1:n])
  eigs = n1 == 1 ? [ A[1,1] ] : eig2(A)

  # swap blocks by repeated rotations
  # invariant: Agnew = (R' * A) * R
  if n == 3
    R = LinAlg.Rotation([g(1, 3)])
    A[:,:] = [ A[3,3] -A[3,2] -A[3,1]; -A[2,3] A[2,2] A[2,1]; -A[1,3] A[1,2] A[1,1]]
  else
    if n1 == 1
      R = LinAlg.Rotation([g(1, 2)])
      A[:,:] = [  A[2,2] -A[2,1]; -A[1,2] A[1,1] ]
    else
      R = g(1, 3) * g(2, 4)
      A[:,:] = [ A[3,3]  A[3,4] -A[3,1] -A[3,2]
                 A[4,3]  A[4,4] -A[4,1] -A[4,2]
                -A[1,3] -A[1,4]  A[1,1]  A[1,2]
                -A[2,3] -A[2,4]  A[2,1]  A[2,2]]
    end
  end

  # perform QR on matrix A
  transform_Hess!(A, 1, n1 + n2, R, eigs, 1)
  reschur!(A, n2)
  A_mul_B!(R', spike)
  spnorm2 = norm(spike[n2+1:n])
  spnorm1, spnorm2, A, spike, R 
end

@inline function reschur!(A::AbstractMatrix, n1::Integer)

  A[n1+1:end,1:n1] = 0
end


"""
Swap two adjacent eigenvalues in diagonal of Schur form.
"""
function swap!(A::AbstractMatrix, ilo::Integer, ihi::Integer, Q::AbstractM, p1::Integer, n1::Integer, p2::Integer, n2::Integer)
  
  n = size(A, 1)
  r1 = p1:p1+n1-1
  r2 = p2:p2+n2-1
  s1 = p2+n2-n1:p2+n2-1
  s2 = p1:p1+n2-1
  A[s1,ilo:n], A[s2,ilo:n] = A[r1,ilo:n], A[r2,ilo:n]
  A[1:ihi,s1], A[1:ihi,s2] = A[1:ihi,r1], A[1:ihi,r2]
  Q[1:ihi,s1], Q[1:ihi,s2] = Q[1:ihi,r1], Q[1:ihi,r2]
  shift = eig2(view(A, r2, r2))
  transform_Hess!(A, p1, p2+n2-1, Q, shift, 2)
  A
end

function eig2(A::AbstractMatrix)
  n = size(A, 1)
  if n == 1
    [ A[1,1] ]
  else
    disc = discriminant(A, 1, 2)
    re = (A[1,1] + A[2,2]) / 2
    sq = sqrt(abs(disc))
    if disc < 0
      [ Complex(re, sq) ]
    else
      [ re - sq; re + sq ]
    end
  end
end

"""
Add missing multiplication to LinAlg.Rotation
"""
function A_mul_Bc!(R::LinAlg.Rotation, G::LinAlg.Givens)
  insert!(R.rotations, 1, G')
  R
end

function A_mul_B!{T<:Number,S<:Number}(A::AbstractMatrix{T}, R::LinAlg.Rotation{S})
  n = length(R.rotations)
  @inbounds for i = 1:n
    A_mul_Bc!(A, R.rotations[n+1-i]')
  end
  A
end

function A_mul_B!(R::LinAlg.Rotation, A::AbstractVecOrMat)
    @inbounds for i = 1:length(R.rotations)
        A_mul_B!(R.rotations[i], A)
    end 
    return A
end

(*){T<:Number,S<:Number}(A::AbstractMatrix{T}, R::LinAlg.Rotation{S}) = A_mul_B!(copy(A), R)

"""
Produce Givens Rotation, which transforms 2x2 matrix to another 2x2 matrix with
either: subdiagonal zero if 2 real eigenvalues.
The eigenvalue of lowest absolute value is placed in lower right position.
or: equal diagonal elements if 2 non-real complex eigenvalues.
The lower left element is <= upper right element absolutely.
"""
function givens1{T<:Number}(A::AbstractMatrix{T}, i1::Integer, i2::Integer)
  a, b, x, d = A[i1,i1], A[i1,i2], A[i2,i1], A[i2,i2]
  btx = b * x
  apd = ( a + d ) / 2
  bpx = ( b + x ) / 2
  da  = ( d - a ) / 2
  disc = da ^ 2 + btx
  if disc >= 0
    root = sqrt( disc )
    root = copysign(root, apd)
    e1 = apd + root
    e2 = ( a * d - btx ) / e1
    ea = a - e2
    ed = d - e2
    G, r = abs(b) + abs(ed) > abs(b) + abs(ea) ? givens(b, ed, i1, i2) : givens(ea, x, i1, i2)
    v = 2( da * G.c - bpx * G.s) * G.s + b
    G, [e1 v; 0 e2]
  else
    bx = ( b - x ) / 2
    root = hypot(da, bpx)
    root = copysign(root, bx)
    w = bx + root
    v = disc / w
    G, r = abs(b) >= abs(x) ? givens(w + x, da, i1, i2) : givens(da, w - b, i1, i2)
    G, [ apd w; v apd]
  end
end

end
