function sylves{T<:Number}(f::AbstractVector{T}, g::AbstractVector{T}, k::Int)
#
#  The k-th sylvester matrix of f(x), f'(x) = g(x)
#
    l = length(f)
    m = l + k - 1
    n = 2 * k + 1
    A = zeros(T, m, n)
	for j = 1:k+1
        A[j:j+l-2,j] = g
    end
    for j = 1:k
        A[j:j+l-1,j+k+1] = f
    end
    A
end
   
function sylves1{T<:Number}(f::AbstractVector{T}, g::AbstractVector{T}, k::Int)
#
#  The k-th sylvester matrix of f(x), f'(x) = g(x)
#
    l = length(f)
    m = l + k - 1
    n = 2 * k + 1
    A = zeros(T, m, n)
	for j = 1:k+1
        A[j:j+l-2,j*2-1] = g
    end
    for j = 1:k
        A[j:j+l-1,j*2] = f
    end
    A
end
   
