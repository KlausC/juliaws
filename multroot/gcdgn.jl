function gcdgn{T<:Number}(p::AbstractVector{T}, q::AbstractVector{T}, g0::AbstractVector{T}, u0::AbstractVector{T}, v0::AbstractVector{T})
#
#  Finds extended GCD of polynomial p and q by Gauss-Newton
#  iteration, such that
#
#         conv(g,u) = p,      conv(g,v) = q
#   
#  Calling syntax:
#    g,u,v,res = gcdgn(p, q, g0, u0, v0)
#      
#        INPUT:  p, q -- polynomial coefficients (hight coeff first)
#                g0, u0, v0 -- initial iterates
#        OUTPUT: g, u, v as described above
#                res
#
	E = real(one(T))
    m = length(g0) - 1; m1 = m+1;   # degree of g (i.e. gcd)
    n = length(u0) - 1; n1 = n+1;   # degree of u
    k = length(v0) - 1; k1 = k+1;   # degree of v
    lp = length(p);                 # length of p
    lq = length(q);                 # length of q = p'
   
    # making all polynomials monic and make a copy
    p  = p / p[1]
    q  = q / q[1]
    g0 = g0 / g0[1]
    u0 = u0 / u0[1]
    v0 = v0 / v0[1]
    
    s = conv(g0, u0) - p
    t = conv(g0, v0) - q
    b = [s[2:lp]; t[2:lq]]   
    x = [g0[2:m1]; u0[2:n1]; v0[2:k1]]
    w = ones(T, length(b))
    for j = 2:lp
	    apj = abs(p[j])
        if apj > E 
            w[j-1] = E / apj
		end   
    end
   
    for j = 2:lq
	    aqj = abs(q[j])
        if aqj > E 
            w[lp-2+j] = E / aqj
        end   
    end
   
    bke0 = norm(abs(b .* w)) 
    bke = bke0 
    #@printf("      g-n %g \n", bke)
	g::AbstractVector{T} = []
	u::AbstractVector{T} = []
	v::AbstractVector{T} = []
   
    j = 1
    while j > 0
       
        A = sylmat(g0, u0, v0)
        d = scalsq(A, b, w)
		y = x
        tau = 1.0
		bkej = bke
		while tau > 1e-4 && bkej > bke * (1-tau*0.25) 
            y = x - d * tau
       
            g = [E; y[1:m]]
	        u = [E; y[m+1:m+n]]
            v = [E; y[m+n+1:m+n+k]]
            s = conv(g, u) - p
		    t = conv(g, v) - q
            b = [s[2:lp]; t[2:lq]]
		    bkej = norm(abs(b .* w)) 
            tau *= 0.5
		end
        #@printf("      g-n %g %g %g \n", bkej, tau, norm(d))
       
        if bkej >= bke ||  bke - bkej < (bke0 - bkej) * 0.05 || tau <= 1e-4
            g = g0
			u = u0
			v = v0
            break
        end

        g0 = g
		u0 = u
		v0 = v
        x = y
		bke = bkej
        j = j + 1;
    end
    g, u, v, bke0, bke
end 

