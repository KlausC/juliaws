"""
Calculation of eigenvalues using Francis Iteration - QR
"""
module Francis

using Base.LinAlg
using LinAlgBigFloat

using LinAlgBigFloat: transform_Hess!, separate!, deflation_criterion, deflation_criterion1, discriminant, seigendiag, givens1, reschur!, finish

import LinAlgBigFloat.AbstractM

"""
Transform a real matrix to Hessenberg form modyfying input matrices.
"""
function hessenberg!{T<:Real}(A::AbstractMatrix{T}, Q::AbstractM)
  transform_Hess!(A, 1, size(A, 1), Q, T[], size(A, 1), 0)
  A, Q
end

"""
Transform a real matrix to Hessenberg form.
"""
function hessenberg{T<:Real}(A::AbstractMatrix{T})
  Q = eye(A)
  hessenberg!(copy(A), Q)
end

"""
Decompose abstract matrix into Schur factors
"""
function schur2(AA::AbstractMatrix)
  n, m = size(AA)
  n == m || error("require square matrix")
  _schur2(AA)
end

function _schur2(AA::AbstractMatrix)
  A = copy(AA)
  Q = eye(A)
  hessenberg!(A, Q)
  separate!(A, 1, size(A, 1), Q, iterate!)
  finish!(A, Q)
  A, Q
end

"""
Find estimations for eigenvalues from lower right of Hessenberg matrix.
"""
function estimations!(A::AbstractMatrix, ilo::Integer, ihi::Integer, Q::AbstractM, iwindow::Integer)
  n, m = size(A)
  n == m || error("square Matrix required")
  ihi = min(n, ihi)
  1 <= iwindow <= ihi || error("submatrix not ok")
  ni = max(ilo-1, ihi - iwindow)
  iwindow = ihi - ni
  if iwindow == 1
    return A
  elseif iwindow == 2
    if discriminant(A, ihi-1, ihi) < 0
      return A
    end
  end
  ra = ni+1:ihi
  rb = ihi+1:n

  # println("before entering smalleig $(ni+1):$ihi")

  A[ra,ra], FZ = smalleig(view(A, ra, ra))
  if ni >= ilo
    A[ra,ni] = FZ' * A[ra,ni]
  end
  A[ra,rb] = FZ' * A[ra,rb]
  A[1:ni,ra] = A[1:ni,ra] * FZ
  Q[:,ra] =  Q[:,ra] * FZ
  A
end

"""
Provide schur matrix and transformation for small matrix
"""
function smalleig(A::AbstractArray)
  # F = schurfact(A)
  # F[:T], F[:Z]
  n = size(A, 1)
  if n == 1
    A, eye(A)
  elseif n == 2
    A, G = givens1(A, 1, 2)
    A, G' * eye(A)
  else
    _schur2(A)
  end
end

"""
Deflate elements in first subdiagonal.
Return new lower and upper indices.
"""
function deflate_subdiagonal!(A::AbstractMatrix, ilo::Integer, ihi::Integer)
  z = zero(eltype(A))
  while (idiff = can_deflate_hi(A, ilo, ihi)) > 0
    ihi -= idiff
    A[ihi+1,ihi] = z
  end
  while (idiff = can_deflate_lo(A, ilo, ihi)) > 0
    ilo += idiff
    A[ilo,ilo-1] = z
  end
  ilo, ihi
end

function can_deflate_hi(A::AbstractMatrix, ilo, ihi)
  if ihi > ilo && deflation_criterion(A[ihi,ihi-1], A[ihi-1,ihi-1], A[ihi,ihi])
    1
  elseif ihi > ilo+1 && deflation_criterion(A[ihi-1,ihi-2], A[ihi-2,ihi-2], A[ihi-1,ihi-1])
    2
  else
    0
  end
end

function can_deflate_lo(A::AbstractMatrix, ilo, ihi)
  if ihi > ilo && deflation_criterion(A[ilo+1,ilo], A[ilo+1,ilo+1], A[ilo,ilo])
    1
  elseif ihi > ilo+1 && deflation_criterion(A[ilo+2,ilo+1], A[ilo+2,ilo+2], A[ilo+1,ilo+1])
    2
  else
    0
  end
end

"""
Early aggressive defaltion step without re-ordering
Return new upper index.
"""
function deflate!(A::AbstractMatrix, ilo::Integer, ihi::Integer, ispike::Integer)
  z = zero(eltype(A))
  n = min(size(A, 1), ihi)
  k = n
  ni = ispike
  if ilo < ispike && deflation_criterion1(abs(A[ilo+1,ilo]), abs(A[ilo,ilo]))
    A[ilo+1,ilo] = z
    ilo += 1
    # println("deflate at top: new ilo = $ilo")
  end
  if ni < ilo
    k = ilo - 1
    # println("deflation ispike = $ni: not required new ihi = $k")
  else
    while k > ni
      akkm = A[k,k-1]
      if akkm == z || k <= ni + 1
        if deflation_criterion1(abs(A[k,ni]), abs(A[k,k]))
          A[k,ni] = z
          k -= 1
        else
          break
        end
      else
        rk = k-1:k
        if deflation_criterion1(norm(view(A, rk, ni)), norm(view(A, rk, rk))) 
          A[k-1,ni] = A[k,ni] = z
          k -= 2
        else
          break
        end
      end
    end
    if k < n
      # println("deflation ispike = $ni: zeroed $(k+1):$n, new ihi = $k")
    end
  end
  ilo, k
end

"""
Repeated aggressive deflation steps re-ordering Eigenvalues
"""
function reorder!{T<:Union{Float64,Float32,Float16}}(A::AbstractMatrix{T}, ilo::Integer, ihi::Integer, Q::AbstractM, iwindow::Integer)

  n = size(A, 1)
  # println("before deflate_sub($ilo,$ihi)")
  ilo2, ihi2 = deflate_subdiagonal!(A, ilo, ihi)
  if ihi2 < ihi
    # println("deflated(A($n $ilo:$ihi)) by $(ihi-ilo-ihi2+ilo2) / $(ihi-ilo+1) to $ilo2:$ihi2")
  end

  ilo, ihi = ilo2, ihi2
  iwindow = min(ihi-ilo+1, iwindow)
  if iwindow <= 0
    return ilo, ihi
  end
  # println("before estimations($ilo,$ihi):"); display(A);
  estimations!(A, ilo, ihi, Q, iwindow)
  ispike = ihi - iwindow
  # println("after estimations($ilo,$ihi): spike at $ispike"); display(A);
  if ispike < ilo
    return ilo, ilo - 1
  end
  # @assert is_hessenberg(A, ispike) "hessenberg after estimations! $ilo:$ihi"
  # @assert is_transform(A, Q) "transform after estimations! $ilo:$ihi"
  k = 0
  loc, hic = ispike+1, ihi
  while loc <= hic && ispike >= ilo
    k += 1
    ilo, ihi2 = deflate!(A, ilo, ihi, ispike)
    # println("after deflate($ilo,$ihi):"); display(A);
    if ihi2 < ihi
      # println("deflated($k A($n $ilo:$ihi)) by $(ihi-ihi2) / $(ihi-ispike) spike $ispike")
      # display(A[ispike+1:ihi,ispike:ihi])
    end
    ihi = ihi2
    hic = min(ihi, hic)
    if loc <= hic
      # loc, hic = loc, loc-1
      loc, hic = swap_sweep!(A, ispike, loc, hic, Q)
      # @assert is_transform(A, Q) "transform after swap_sweep $ilo:$ihi"
    end
  end
  if ihi >= ilo
    ev = seigendiag(A, ispike+1, ihi)
    transform_Hess!(A, ilo, ihi, Q, zeros(T, 0), ihi, 0) # remove spike from A
    # @assert is_hessenberg(A) "removed spike at $ispike - $ilo:$ihi"
    for i = 1:min(max(1,(ihi-ilo+1)÷3),3)
      for k = length(ev):-1:1
        transform_Hess!(A, ilo, ihi, Q, [ev[k]], 1, 0) # insert eigenvalue estimations
      end
    end
    transform_Hess!(A, ilo, ihi, Q, zeros(T, 0), ihi, 0) # remove bulges from A
    # @assert is_hessenberg(A) "removed bulges $ilo:$ihi"
    # @assert is_transform(A, Q) "transform after bulges chase $ilo:$ihi"
  end
  ilo, ihi
end

"""
Iterate reorder! till convergence
"""
function iterate!(A::AbstractMatrix, ilo::Integer, ihi::Integer, Q::AbstractM)

  while ilo <= ihi
    iwindow = window_size(A, ilo, ihi)
    # println("reorder!(A, $ilo, $ihi, Q, $iwindow)")
    ilo, ihi = reorder!(A, ilo, ihi, Q, iwindow)
    # @assert is_hessenberg(A) "is_hessenberg $ilo:$ihi"
    # @assert is_transform(A, Q) "is_transform $ilo:$ihi"
  end
end

"""
Find appropriate window size. Must be >= 2 and < ihi - ilo + 1.
"""
function window_size(A, ilo, ihi)
  if ihi - ilo <= 1
    return 2
  end
  wlo, whi = window_size_heuristic(A, ilo, ihi)
  wlo = max(wlo, 2)
  whi = min(whi, ihi - ilo)
  while wlo >= whi
    whi += 1
  end
  whi = min(whi, ihi-ilo)
  # find small subdiagonal element in interval [jlo, jhi]
  best = Inf
  wbest = 2
  for w = wlo:whi
    bestw = abs(A[ihi-w+1,ihi-w])
    # print("w = $w A[$(ihi-w+1);$(ihi-w)] = $bestw")
    if bestw < best
      wbest = w
      best = bestw
      # print(" ***")
    end
    # println()
  end
  wbest
end

"""
Lower and upper limits for window size.
"""
function window_size_heuristic(A, ilo, ihi)
  k = min(256, ( ihi - ilo + 1 ) ÷ 3 )
  (k * 9) ÷ 10, (k * 11 + 9) ÷ 10
  # 2, (k * 11 + 9) ÷ 10
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
Perform swap on copy of small submatrix and check if spike vector is improved by the swap.
Return true if last component of changed spike vector is smaller than before.
and return range of transformation, changed small matrix, spike vector, transformation matrix.
"""
function testswap(A::AbstractMatrix, spike::AbstractVector, ilo::Integer, iup::Integer)
  i1, n1, n2 = index_sizes(A, ilo, iup)
  if i1 >= ilo
    r = i1:iup
    sp1, sp2, A2, spike2, R = swap_small!(A[r,r], spike[r], n1, n2)
    if sp2 < sp1
      true, r, A2, spike2, R
    else
      return false, r, A2, spike2, R  
    end
  else
    false, i1:i1-1, 0, 0, 0 
  end
end

"""
Starting from end, repeatedly swap adjacent eigenvalue blocks until begin of diagonal.
"""
function swap_sweep!(A::AbstractMatrix, ispike::Integer, ilo::Integer, ihi::Integer, Q::AbstractM)
  n = size(A, 2)
  i = ihi
  ilo = max(ilo, ispike + 1)
  hic = 0
  loc = n + 1
  # println("swap_sweep($ispike, $ilo:$ihi)")
  while i >= ilo
    res, ra, A2, spike2, R = testswap(A, A[:,ispike], ilo, i)
    # println("testswap($ilo,$i): $res $ra $spike2")
    if length(ra) > 0 || res 
      rb = last(ra)+1:n
      rn = 1:first(ra)-1
      A[ra,ispike] = spike2
      A[ra,ra] = A2
      A[ra,rb] = R' * A[ra,rb]
      A[rn,ra] = A[rn,ra] * R
      Q[:,ra] = Q[:,ra] * R
      loc = first(ra)
      if hic == 0
        hic = last(ra)
        if hic < ihi
          hic += 1
          if hic < ihi && A[hic+1,hic] != 0
            hic += 1
          end
        end
      end
    end
    i = first(ra)
    if i >= ilo
      if A[i+1,i] != 0
        i += 1
      end
    end
  end
  if loc > ispike + 1
    loc -= 1
    if loc > ispike + 1 && A[loc,loc-1] != 0
      loc -= 1
    end
  end
  # loc, hic, A, Q
  loc, loc-1, A, Q
end


"""
A is a n1+n2:n1+n2 Schur matrix, spike a n1+n2 vector. 1 <= n1, n2 <= 2.
Orthogonally transform A to a Schur matrix with the diagonal eigenvalue blocks swapped.
"""
function swap_small!{T<:Real}(A::AbstractMatrix{T}, spike::AbstractVector{T}, n1::Integer, n2::Integer)

  ze, on = zero(T), one(T)
  g(i1::Integer, i2::Integer) = givens(ze, on, i1, i2)[1]
  n = n1 + n2 
  # source ranges
  spnorm1 = norm(spike[n1+1:n])
  eigs = n1 == 1 ? [ A[1,1] ] : eig2(A)

  # swap blocks by repeated rotations
  # invariant: Anew = (R' * A) * R
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
  transform_Hess!(A, 1, n1 + n2, R, eigs, 1, 0)
  reschur!(A, n2)
  A_mul_B!(R', spike)
  spnorm2 = norm(spike[n2+1:n])
  spnorm1, spnorm2, A, spike, R 
end

"""
Test if square matrix A is transformed to A by Q. (B * Q - Q * A) is small.
"""
function is_transform(B::AbstractMatrix, A::AbstractMatrix, Q)
  norm(Q * A - B * Q, 1) <= eps(eltype(A)) * size(A, 1) * 1000
end

end

