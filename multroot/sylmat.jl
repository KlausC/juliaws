function sylmat{T<:Number}(g::AbstractVector{T}, u::AbstractVector, v::AbstractVector{T})
#
#  sylmat generates the Sylvester resultant matrix
#    of polynomial g, u, v
#
    mg = length(g) - 1 
    mu = length(u) - 1 
    mv = length(v) - 1 
    np = mg + mu
    nq = mg + mv
   
    m = mg + mu + mv
    n = np + nq
      
    A = zeros(T, n, m)  
   
    for j = 1:mg
       A[j:j+mu,j] = u
       A[np+j:np+j+mv,j] = v
    end
   
    for j = 1:mu
       A[j:j+mg,j+mg] = g
    end
   
    for j = 1:mv
       A[np+j:np+j+mg, j+mg+mu] = g
    end
    A
end
