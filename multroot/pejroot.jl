function pejroot{T<:Number,S<:Number,U<:AbstractFloat}(f::AbstractVector{T}, z0::PolyZeros{S},
					  noi::Int = 10, tol::U = 1e-8, style::Int = 2, prtsty::Int = 0)
#
# PEJROOT calculates (multiple) roots of polynomial using an iterative
# method developed by Z. Zeng. For details, you can request the paper
#
#     "Conputing Multiple roots of polynomials: pejorative condition and accurate
#     computation", by Zhonggang Zeng 
#
# This code is dated Dec. 4, 2002. Please sent bug/failure report to
#           
#             Prof. Zhonggang Zeng
#             Department of Mathematics
#             Northeastern Illinois University
#             Chicago, IL 60625
#
#             email: zzeng@neiu.edu
#
# This code is freely released for research exchange only. The author 
# is not responsible for any demage caused by using this code.
#
# Calling syntax: The simplest way to call is
#  
#   » pejroot(f,y,m)
#
# where the input items are
#          f --- (row vector) the target polynomial
#          y --- (row vector) the initial iterate (of the roots)
#          m --- (row vector) the multiplicity structure
#
# For more advanced usage:
#
#	» [z, e, c] = pejroot(f, y, l, noi, tol, style, prtsty)
#
# The output 
#          z --- (matrix) distinct roots (1st column) and corresponding
#                         multiplicities (2nd column)
#          e --- (scalar) the backward error of the roots
#          c --- (scalar) the pejorative condition number
#
#
#  Additional input parameters
#  
#        noi --- (integer)    the number of iteration allowed
#                               (default: 10)
#        tol --- (numeric)    the error tolerance
#                               (default: 1.0d-8)
#      style --- (integer)    backward error minimization style:
#                               1 : overall (absolute)
#                               2 : coefficientwise (relative,default)
#     prtsty --- (integer)    intermediate results showing style
#                               0: minimal (default)
#                               1: plus intermediate backward error and 
#                                     root corrections
#                               2: plus intemediate root approximations
#  Example:
#
# » f = poly([ones(1,20) 2*ones(1,10) 3*ones(1,5)]); # construct the test 
#                					# polynomial with multiple roots
# » y = [0.95,2.05,2.95];           # prepare initial iterate
# » m = [20, 10, 5];                # prepare the multiplicity structure
# » pejroot(f,y,m)                  # running the program
# ans =
#
#   3.00000000000000   5.00000000000000
#   2.00000000000000  10.00000000000000
#   1.00000000000000  20.00000000000000
#
#
    CT = coeffstype(z0.z)
    if ! issubtype(CT, T) throw(ArgumentError("coefficient type $T not sufficient for $CT")) end
    E = one(U)
    Z = zero(U)
    m = length(z0.mult)	# number of variables
    n = sum(z0.mult)		# number of equations
    f = f / f[1]		# make the polynomial monic
   
    if length(f) != n+1   # exit if the input is wrong
        throw(ArgumentError("Input error"))
    end
    
    # sort initial values to enhance accuracy
    # It is interesting, sorting really improves accuracy
    jj = sortperm(abs(z0.z), rev=true)	# greatest first
    y0 = z0.z[jj]
    l = z0.mult[jj]

    job = 0						# initialize job
    y = conj(y0)				# pass the initial iterate
    h = conj(f[2:n+1])			# make the RHS
    delta = zeros(U, 1, noi)	# space for sizes of the correction
    bcker = zeros(U, 1, noi)	# space for backward errors

    w = ones(T, n)				# set weight
    if style == 2
       for j = 1:n
           afj = abs(f[j+1])
           if afj >= E 
               w[j] = E / afj
           end
       end
    end

    if prtsty == 1 
       println("     step   bcker         delta")
    end

    Df = zeros(T, n ,m)			# open the space for Df
    bkerr = Z 
    pjcnd = Z
    k = 0

    for k = 1:noi

        # evaluate the coefficient operator and its Jacobian
	println("before Zeros")	
        z = PolyZeros(y, l)
	println("after Zeros")	
        b = polyG(z) - h
	println("before polyDG")	
        Df = polyDG(z)
	println("after polyDG")	
      
        if style == 2   # scale
            b = b .* w
            for j = 1:n
                Df[j,1:end] = Df[j,1:end] * w[j]
            end
         end

         d = Df \ b							# least squares 
         delta[k] = norm(d,2)
         bcker[k] = norm(b,Inf)

         if delta[k] < tol
            job = 1
         end     							# convergence criterion 1
         if k > 1 
             if delta[k] > delta[k-1] && bcker[k] > bcker[k-1]
                 if job == 1
                     bkerr = bcker[k-1]
                     break
                 end
             elseif delta[k] < delta[k-1]	# criterion 2
                 if delta[k]^2 / (delta[k-1] - delta[k]) < tol 
                     job = 1 
                 end
             end
         end
      
         y = y - d		# correct the roots

         if prtsty >= 1		# intermediate show style 1
             @printf(" %2.0f    %10.2e    %10.2e", k, bcker[k],delta[k])
         elseif prtsty == 2	# intermediate show style 2
             println(y)
         end
  
          bkerr = bcker[k]
      end # k-loop

      if job == 1
          s = svd(Df)
          pjcnd = 1 / s[2][m]  # get pej. cond. number
          jj = sortperm(l)
          y = y[jj]
          l = l[jj]   # sort by multiplicities
      end

      PolyZeros(y, l), bkerr, pjcnd, job, k
end
