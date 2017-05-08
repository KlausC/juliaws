function scalsq{T<:Number}(A::AbstractArray{T,2}, b::AbstractVector{T}, w = nothing)
#
#  Solving the scaled least squares problem 
#          W(Ax - b) = 0
#      with iterative refinement. 
#
#  input  A --- the matrix
#         b --- the right-hand side vector
#         w --- the scaling vector (diagonal entries of W)
#
#      
    m, n = size(A)
	E = real(one(T))
   
    if w == nothing
        w = ones(T, m)
        for j = 1:m
            abj = abs(b[j])
            if abj > E
                w[j] = E / abj
            end
        end
    end
   
    for j = 1:m
        A[j,1:end] = A[j,1:end] * w[j]
    end

    b = b .* w
   
    Q, S = qr(A)

    d = Q' * b
    x = backsub(S, d)
   
    # one step refinement
    bb = [b; zeros(T, n)]
	B = [eye(T, m) A; A' zeros(T, n, n)]
    r = b - A * x
   
    #for j = 1:3
    rr = bb - B * [r; x]; #disp(norm(rr));
   
    s = Q' * rr[1:m]
    c = forsub(S', rr[m+1:m+n])
	c2 = backsub(S, s[1:n] - c)
    #c1 = Q * [c; s[n+1:m]]
    #r = r + c1
    x = x + c2
    #%end;
    x
end
