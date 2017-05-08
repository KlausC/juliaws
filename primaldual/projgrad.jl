
function penalty_functions()
#=
    g(x) = sum(log(1.0 ./ max(x, 0.0)))
    Dg(x) = - 1.0 ./ max(x, 0.0)'
    D2g(x) = (1.0 ./ max(x, 0.0) .^ 2)
    g(x) = sum(1.0 ./ max(x, 0.0))
    Dg(x) = (- 1.0 ./ max(x, 0.0) .^ 2 )'
    D2g(x) = (2.0 ./ max(x, 0.0) .^ 3)
=#
    k = 1
    s(x) = sum(max(x, 0.0) .^ -k)
    r(y::Float64) = y ^ (1/k)
    r1k(y::Float64) = r(y) / y
    r2k2(y::Float64) = -(k - 1) * r1k(y) / y
    g(x) = r(s(x))
    Dg(x) = - r1k(s(x)) * ( max(x, 0.0) .^ (-k-1))'
    function D2g(x)
      sx = s(x)
      [ r1k(sx) * (k + 1) * (max(x, 0.0) .^ (-k-2))   sqrt(-r2k2(sx)) * (max(x, 0.0) .^ (-k-1)) ]
      # [ d a ] represents matrix Diagonal(d) - a * a'
    end
	g, Dg, D2g
end # functions

function dirfuncs(x::Vector)
   
    A = Dh(x)
    m, n = size(A)
    F = lufact(A)   # selection of base variable is left to magic of lufact
    LL, UU, p, q, Rs = F[:(:)]
    pin = sortperm(p)
    qin = sortperm(q)
    q1 = q[1:m]
    q2 = q[m+1:n]

    L = LowerTriangular(LL)
    U1 = UpperTriangular(UU[:,1:m])
    U2 = UU[:, m+1:n]
    
    # Tangential mapping.  ( Dh * T = 0)
    T(y)  =  [ - U1 \ (U2 * y); y][qin,:]

	# Transposed of T operation on row vectors			 
    Tl(r) = r[:,q2] - (r[:,q1] / U1) * U2 

	# Map residual value of equations to base variables (q1)
    B(h) = U1 \ (L \ ((Rs .* h)[p,:]))

	# Restauration direction (only q1-Variables != 0)
    function d(h)
	    dd = spzeros(n, size(h,2))   
		dd[q1,:] = B(h)
		dd[:,1]
    end

	# Lagrange parameters with respect to linear equations A
    lambda(df) = -(df[:,q1] / U1 / L)[:,pin] .* Rs'

	reg = 1.0	# added to diagonal of D2g to avoid quasi singularity
	stol = 1e-4 # reduction factor for residuum for pcg method
	mit = min(200, (n-m) ÷ 50) # maximal iteration count for pcg method

	# from here L denotes the second derivative of function (only comments)
	dg = D2g(x)   # second derivative of taget function plus regularizer
	d1 = dg[q1,1] + reg # diagonal part - base variables
	d2 = dg[q2,1] + reg #               
	if size(dg, 2) == 2
      e1 = dg[q1,2] # rank 1 modification - base variables 
	  e2 = dg[q2,2] #      
	  eU = e2 - U2' * ( U1' \ e1) # -eU * eU' is the rank1 modification of TLT 
	  cd = d2 + ((U1 \ U2)' .^ 2) * d1  # diagonal of first part of A
	  da = eU ./ cd
	  det = 1.0 - vecdot(eU, da)
	  # Preconditioner: cheap approximation to inverse of TLT
	  # PRE(y) =  y ./ (cd - eU .^2)
	  PRE(y) = (vecdot(da, y) / det) * da + y ./ cd # preconditioner is (diagonal - rank1) ^-1
	  # Apply T * L * T to vector
	  # A = TLT is symmetric matrix  (D2 + U2' U1'^-1 D1 U1^-1 U2) - (e2 - U2' U1'^-1 e1)(...)'
	  function TLT(y)
	    Uy = U1 \ (U2 * y)
	    U2' * (U1' \ (d1 .* Uy)) + d2 .* y - eU * (vecdot(e2,y) - vecdot(e1,Uy))
	  end
    else
	  eU = 0.0
	  e1 = d1
	  e2 = d2
	  PRE(y) = y ./ cd # preconditioner is diagonal
	  function TLT(y)
	    Uy = U1 \ (U2 * y)
	    U2' * (U1' \ (d1 .* Uy)) + d2 .* y
	  end
	end

    # approximately apply inverse of TLT to right hand side vector using pcg-method
    function TLTm1(rhs)
	    x0 = PRE(rhs)
		nr, kmax, x, r = pconjgrad(x0, TLT, PRE, rhs, mit, stol)
		nr <= stol ? x : x0
	end

	# Apply T * (T' L T)⁻1 * T' to df
	function TTLTT(df)
	    rhs = vec(Tl(df))
		y = TLTm1(rhs)
		z = T(y)
		z[:,1]
	end

	B, T, Tl, d, lambda, TLT, TLTm1, TTLTT, PRE
end
