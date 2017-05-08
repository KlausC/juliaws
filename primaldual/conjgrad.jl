
#=
In most cases, oner matrix M has to be symmetric positive-definite and fixed, i.e., cannot change from iteration to iteration. If any of these assumptions on the preconditioner is violated, the behavior of the preconditioned conjugate gradient method may become unpredictable.preconditioning is necessary to ensure fast convergence of
the conjugate gradient method.
The preconditioned conjugate gradient method takes the following form:

    r0 := b − A x0
    z0 := inv(M) r0
    p0 := z0
    k := 0
    repeat

        αk := rk T zk / pk T A pk
        xk+1 := xk + αk pk
        rk+1 := rk − αk A pk
        if rk+1 is sufficiently small then exit loop end if
        zk+1 := inv(M) rk+1
        βk := zk+1 T rk+1 / zk T rk
        pk+1 := zk+1 + βk pk
        k := k + 1
    end repeat
    The result is xk+1
=#
#=
The above formulation is equivalent to applying the conjugate gradient method without
preconditioning to the system

    inv(E) A inv(E) T x^ = inv(E) b
where
    E E T = M
    x^ = E T x

The preconditioner matrix M has to be symmetric positive-definite and fixed, i.e.,
cannot change from iteration to iteration. If any of these assumptions on the preconditioner
is violated, the behavior of the preconditioned conjugate gradient method may become
unpredictable.

=#

function pconjgrad{T <: AbstractFloat}(x0::AbstractVector{T}, A, C, b::AbstractVector{T}, mit::Int, stol::T)

  fA = typeof(A) <: Function ? A : x -> A * x  
  fC = typeof(C) <: Function ? C : x -> C * x

  nb = norm(b)
  xkk = x0
  rkk = b - fA(xkk)
  zkk = fC(rkk)
  pkk = zkk
  nrkk = norm(rkk) / nb
  k = 0
  println("k: ", k, " rkk: ", nrkk)

  k = 1
  while k <= mit && nrkk > stol
    xk, pk, rk, zk, nrk = xkk, pkk, rkk, zkk, nrkk
    Apk = fA(pk)
    ak = vecdot(rk, zk) / vecdot(pk, Apk)
    xkk = xk + ak * pk
    rkk = rk - ak * Apk
    nrkk = norm(rkk) / nb
    delta = -ak * (pk' * rkk * 2 + pk' * Apk)
    println("k: ", k, " rkk: ", nrkk, " ðtf: ", delta, " rkk/rk: ", nrkk / nrk)
    nrkk > stol || continue

    zkk = fC(rkk)
    #bk = vecdot(zkk, rkk) / vecdot(zk, rk)
    bk = vecdot(zkk, rkk - rk) / vecdot(zk, rk)
    pkk = zkk + bk * pk
    k = k + 1
  end
  
  nrkk, k, xkk, rkk

end
