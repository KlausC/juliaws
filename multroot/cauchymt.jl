function cauchymt{T}(f::AbstractVector{T}, k::Int)
#
# generate a k-column Cauchy matrix of polynomial f
#
    n = length(f)
    m = n + k - 1
    A = zeros(T, m, k)
   
    for j = 1:k
        A[j:j+n-1,j] = f
    end
    A
end
   
