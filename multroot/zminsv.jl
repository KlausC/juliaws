function zminsv{S<:Real,U<:AbstractFloat}(A::AbstractArray{S,2}, tol::U)
#
#  zminsv calculates the smallest singular value
#  of matrix A and the associated right singular vector
#  via an implicit inverse iteration in the form of the
#  Gauss-Newton iteration
#
#    input  A --- the matrix
#           tol --- the error tolerence
#
#   output  s --- the smallest singular value
#           x --- the associated right singular vector
#
    E = one(S)
    m, n = size(A)           	# get the dimensions of A
    if m < n-1 throw(ArgumentError("zminsv only if m >= n-1 but $m < $n-1")) end
   
    Q, R = qr(A)             	# QR decomp. of A, maybe input
    if m == n - 1				# Corner case when s == 0 is guranteed 
        y = normalize([backsub(R[1:m,1:m], R[1:m,end]); -E])
        ss = zero(S)
	    # println("zminsv (m=$m, n=$n) returns ss= $ss y = $y")
        return ss, y
    end

    cr = 0
    tau = norm(A, Inf)     			# get the magnitude of rows as scaler
	x = normalize(map(S, rand(n)) * 2 - 1)		# random initial vector (row)
	r = R * x
    xn = E
	y = x

    for j = 1:5
		D = [ x' * tau * 2; R]		# Least squares problem min |D z - b|²
		b = [ (xn * xn - E) * tau; r ]
        QQ, RR = hessqr(D) 			# Hessenberg QR decomp. of stacked matrix
        z = hqrt(QQ, b)        		# same QQ on b
   
        z = backsub(RR[1:n,1:n], z[1:n])	# backward substitution to solve RR z = QQ' b
        y = x - z					# new iteration
		xn = norm(y)
        r = R * y
        ss = norm(r) / xn
        crk = norm(y-x)
		# println("ss = $ss crk = $crk")
      
        if (j == 1 && crk < tol ) || (j > 1 &&  crk < cr && crk^2 / (cr - crk) < tol) 
            break
        end
        x = y
        cr = crk
		# println("hessenberg least square iteration $j: smin = $ss")
    end
	x = y / xn
	# println("zminsv (m=$m, n=$n) returns ss= $ss x = $x")
    ss, x
end


if ! isdefined(:normalize)
function normalize(v::AbstractVector, p::Real = 2)
    invnrm = inv(norm(v, p))
	vv = copy(v)
	scale!(vv, invnrm)
	vv
end
end


