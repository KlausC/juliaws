function backsub{T<:Number}(A::AbstractArray{T,2}, b::AbstractVector{T})
#
# assume A is an upper triangular matrix
# backsub performs backward substitution to solve Ax=b
#
    n = size(A,1)
    N = zero(T)
    s = maximum(abs(diag(A))) * eps(T)
    x = zeros(T, n)
   	if A[n,n] == N A[n,n] = s end
    x[n] = b[n] / A[n,n]
    for k = (n-1):-1:1
        akk = A[k,k]
        if akk == N akk = s; A[k,k] = s end
        x[k] = (b[k] - vecdot(A[k,k+1:n], x[k+1:n])) / akk 
    end
    x
end

