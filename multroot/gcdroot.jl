function gcdroot{T<:Number,S<:Real}(p::AbstractVector{T}, tol::S = 1e-10)
#  
#  GCDROOT calculates the multiplicity structure and initial root 
#  approximation of a real polynomial p,
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
#  INPUT :  p = polynomial coefficients. 
#           tol = zero remainder threshold; 
#                 default = 10^{-10}, set internally
#                 if no value is specified.
#
#  OUTPUT : z = roots of polynomial p
#           l = corresponding multiplicities
#           bkerr = backward error
#
#  CALL :   It is best to call pzero with only one argument as   
#               >> z = gcdroot(p). 
#
    gamma = S(100)		# residual growth factor, default = 100
    delta = S(100.0)	# threshold growth factor, default = 100
    thresh = tol*100	# zero singular threshold, default = 100*tol
    drop = S(5.0e-5)	# try G-N if sigma(k) < drop * sigma(k-1)
    E = one(T); RE = real(E); Z = zero(T); RZ = real(Z)

    # make the polynomial monic
    if p[1] == Z
        p = p[find(p):end]
    end
    p = p / p[1]
    n = length(p) - 1
    q = deriv(p) / n #  q(x) = p'(x) / degree(p)

    f = copy(p); g = copy(q)		# back up the polynomials 
    nf = maxabs(f)	 				# the largest coefficient

    mx = n; wtol = tol; s0 = RZ; s = RZ; wthrh = thresh

    k = n;            # the degree of working polynomial
    while k >= 1 
        if k == 1     # the polynomial is linear, GCD = 1
            h = E; u = f; v = E; m = 1
        else
            for m = 1:k
                A = sylves1(f, g, m)
                scalerows!(A)
                s0 = s
            
                s, x = zminsv(A, tol)
                if k <= n
					analysis_sv(k, m, A, x, s)
                end
                # @printf("s. value %g,%g,%g\n", m, k, s)
                if s < wthrh * nf || m == mx || s < drop * s0
                    h, u, v, res0, res, sm = gcd_refinement(x, m, f, g, A)
                    if res < wtol || m == mx
                        wtol = max(wtol, res * gamma)	# increase tolerance by factor 
                        wthrh = max(wthrh, (s / nf) * delta)	# increase threshold by factor
                        break # for m
                    end
                end
            end # for m
        end
        # println("after m-loop: m = $m mx=$mx k=$k n=$n s=$s thr=$(wthrh*nf) drop=$(drop*s0) ")
        if k == n			# the root values of u contain all roots of f 
            z = roots(u)	# u has only simple roots
            l = ones(Int64, m) 
            if m == 1
                return finish(z, l, p)
            end
            #horner_analysis(z, f)
        else
            t = roots(u)	# u has only simple roots
            jj = 0
            for j = 1:m
                tj = t[j]
                _, jj = findmin(abs(z - tj))	# find root closest to tj
                ljp = l[jj] + 1
                l[jj] = ljp
				# z[jj] += (tj - z[jj]) / ljp		# store mean value with weights l[jj] and 1
            end
            if m == 1
                l[jj] = l[jj] + k - 1
            end
        end
        if m > 1
            k = k - m
        else
            k = 0
        end
        if k <= 0
            return finish(z, l, p)
        end
        f = h
        g = deriv(f) / k
        nf = norm(f)
        mx = m
    end # k-while loop
end

function finish(z, l, p)
    f = polymult(z, l)
    w = ones(eltype(p), length(p))
    for j = 1:length(p)
        apj = abs(p[j])
        if apj > 1
            w[j] = 1 / apj
        end
    end
    # println("f = $(size(f)) p = $(size(p)) w = $(size(w))")
    bkerr = maxabs((f - p) .* w)
   
    PolyZeros(z, l), bkerr
end

function gcd_refinement(x, m, f, g, A)
    u0 = x[1:2:end] / x[1]
    v0 = x[2:2:end] / x[2]	# u(x) * g(x) == v(x) * f(x)
    #println("refinement u0: $(norm(conv(u0,g) - conv(v0,f))/norm(f))")
    B = cauchymt(u0, length(f) - length(u0) + 1)
    g0 = scalsq(B, f)	# scaled least squares solver g0 = div(f, u0)
    g0 = g0 / g0[1]		# first approximation of gcd(p, p')
    #printresidue("g0", g0, u0, v0, f, g)
    h, u, v, res0, res = gcdgn(f, g, g0, u0, v0)	# refinement of g, u, v by G-N
    #println("refinement u: $(norm(conv(u,g) - conv(v,f))/norm(f))")
    #printresidue("h ", h, u, v, f, g)
    sm = norm(A * [u; v]) / norm([u;v])
    h, u, v, res0, res, sm
end

function printresidue(text, g0, u0, v0, f, g)
    uuu = conv(g0, u0) - f	# should be zero
    vvv = conv(g0, v0) - g	# should be zero

    # println("$text p(x)  $f")
    # println("$text p'(x)  $g")
    println("$text u(x)  $u0")
    println("$text v(x)  $v0")
    # println("$text h(x)  $g0")
    # println("$text Δp(x)  $uuu")
    # println("$text Δp'(x) $vvv")
end

function analysis_sv(k, m, A, x, s) 

    #println("svcheck: $k-$m norm(Ax) = $(norm(A*x)) s = $s")


end

function scalerows!(A::Array{Float64,2})
	s = map(i->2.0^-exponent(norm(A[i,:])), 1:size(A,1))
    scale!(s, A)
	s
end

