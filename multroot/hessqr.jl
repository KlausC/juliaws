function hessqr{T<:Number}(A::AbstractArray{T,2})
#
#  hessqr performs QR decomposition on an upper
#      Hessenberg matrix
#
#  assume A is a Hessenberg matrix
#  output  B -- upper triangular
#          t -- rotation used
#
   m, n = size(A)
   if m < n throw(ArgumentError("m ($m) must be >= n ($n)")) end
   
   R = copy(A)
   Q = zeros(T, 2, min(n,m-1))
   
   for j = 1:n
       if j < m
           rjj = R[j,j]; rjpj = R[j+1,j]
           d = hypot(abs(rjj), abs(rjpj))
           if d != 0
               c = conj(rjj) / d; s = conj(rjpj) / d
               Tr = [c s; -conj(s) conj(c)]
               R[j:j+1,j:n] = Tr * R[j:j+1,j:n]
               Q[1,j] = c; Q[2,j] = s
           else
               Q[1,j] = 1; Q[2,j] = 0
           end
       end
   end
   Q, R
end
