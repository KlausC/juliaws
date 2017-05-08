
function hqrt(t, b)
#
#  transformation using the t from hessqr
#
   z = copy(b)
   k = min(size(t,2), size(z,1)-1)
   for j = 1:k
       T = [t[1,j] t[2,j]; -conj(t[2,j]) conj(t[1,j])]
       z[j:j+1,1:end] = T * z[j:j+1,1:end]
   end
   z
end

function hqrt1{T}(t::AbstractArray{T,2}, n::Int=1)
#
#  transformation using the t from hessqr on unit vector e(n)
#
   k = size(t,2)
   z = zeros(T, k+1)
   z[n] = one(T)
   for j = n:k
	   zj = z[j]
	   z[j] *= t[1,j]
	   z[j+1] = t[2,j] * zj
   end
   z
end
