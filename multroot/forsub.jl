function forsub{T<:Number}(A::AbstractArray{T,2}, b::AbstractVector{T})
#
# assume A is an lower triangular matrix
# forsub performs forward substitution to solve Ax=b
#
    n = size(A,1)
    x = zeros(T, n)
   
    x[1] = b[1] / A[1,1]
    for k = 2:n
        x[k] = (b[k] - vecdot(A[k,1:k-1], x[1:k-1])) / A[k,k]
    end
    x
end
