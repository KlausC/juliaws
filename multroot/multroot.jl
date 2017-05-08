function multroot{T<:Number,S<:Real}( p::AbstractVector{T}, tol::S = 1e-10)
#  
#  Finds all roots of a real or complex polynomial p,
# 
#        p(x) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_n+1,
#
#  given by an (n+1)-row vector p = (a_1, ..., a_n+1) 
#
#  METHOD : for details, contact the author 
#
#             Zhonggang Zeng
#             Department of Mathematics
#             Northeastern Illinois University
#             Chicago, IL 60625
#          
#             email: zzeng@neiu.edu
#
# This code is freely released for research exchange only. The author 
# is not responsible for any demage caused by using this code.
#
#  INPUT :  p = polynomial coefficients. 
#           
#  OUTPUT : z = distinct roots of polynomial p
#           l = corresponding multiplicities
#  CALL :   
#               >> z = multroot(p). 
#
    #
    # clear leading/trailing zeros
    #
    n = length(p)
	Z = zero(T)
	E = one(T)
    if p[1] == Z || p[n] == Z 
       jj = find(p) 
       j1 = minimum(jj)
	   j2 = maximum(jj)
       q = p[j1:j2]
    else
        j1 = 1; j2 = n;
        q = p
    end
    q = q / q[1]
    #
    # scaling
    #
    m = length(q)-1
    c = E / (abs(q[m+1])) ^ (1/m)
    # println("c = $c, m = $m, q = $q")
   	q = q .* (c .^ (0:m))
println("before gcdroot")
	z0, bke = gcdroot(q, tol)
   
    if bke < S(1.0e-2)
println("before pjeroot")
        z1, bkerr, pjcnd, job = pejroot(q, z0)
        z = z1.z
        l = z1.mult
        if j2 < n
            z = [z; Z]
			l = [l; n-j2]
        end
        if job == 1
            #
            # show off results
            #
            z = z / c
			nz = length(z)
            @printf("\n");
            @printf("    !!!THE COMPUTATION IS SUCCESSFUL!!!\n")
            @printf("\n")
            @printf("THE PEJORATIVE CONDITION NUMBER:       %g \n", pjcnd)
            @printf("THE BACKWARD ERROR:                    %6.2e \n", bkerr)
            @printf("THE ESTIMATED FORWARD ROOT ERROR:      %6.2e \n", 2 * bkerr * pjcnd)
            @printf("\n");
            if norm(imag(z)) == Z 
                @printf("        computed roots         multiplicities\n")
                @printf("\n")
				for j = 1:nz @printf("%25.15f \t \t \t %3g \n", z[j], l[j]) end
            else
                @printf("        computed roots ")
                @printf("   \t\t\t\t\t\t     multiplicities\n");
                @printf("\n");
                for j = 1:nz
				    @printf("%22.15f + %22.15f i \t \t %3g \n", real(z[j]), imag(z[j]), l[j])
				end
            end
        else
            z = roots(pp)	# fallback to standard rootfinder
			l = ones(n - 1)
        end
    else
        z = roots(pp)		# fallback to standard rootfinder
		l = ones(T, n - 1)
		job = 0
		bkerr = S(Inf)
		pjcond = S(Inf)
    end
    z1, bkerr, pjcnd, job
end

