

# given an array a[1:n] of real numbers, find three subarrays a[1:k1], a[k1+1;k2], a[k2+1;n]
# in a way that minimizes the difference between a and constants on subarrays.
#

type Sums
    a::AbstractArray{Float64,1}
	aa::AbstractArray{Float64,1}
	aaa::AbstractArray{Float64,1}
	perm::AbstractArray{Int64,1}
	optind::AbstractArray{Int64,1}
	optval::Float64
	optx::AbstractArray{Float64,1}
	function Sums(a)
	    n = length(a)
        p = sortperm(a, rev=true)
        b = a[p]
		aa = cumsum(b)
	    new(b, aa, cumsum(aa), p)
    end
end

function sums(s::Sums, k::AbstractArray{Int64,1})
    m = length(k)
	n = maximum(k)
	r = 0.0
	kk = max(1, sort(map(Int, k)))
	for i = 2:m
        k2 = min(max(kk[i],1),n)
		k1 = min(max(kk[i-1],1),n)
		x = s.aa[k1]
		y = s.aa[k2]
		r += (y + x) * (k2 - k1)
    end
	r += s.aa[kk[1]] + s.aa[kk[end]]
	r * 0.5
end
    
function deviation(s::Sums, k)
    n = minimum(k)
	m = maximum(k)
	x = (n > 1) ? s.aaa[n-1] : 0.0
	int = s.aaa[m] - x
	1.0 - sums(s, k) / int
end	

function alldeviations(s::Sums)
    n = length(s.a)
	for k1 = 1:n
        for k2 = k1:n
            s.all[k1,k2] = deviation(s, [1,k1,k2,n])
        end
    end
end

function gradstep(x, f, n)
    d = length(x)
	g = zeros(d)
    fx = f(x)
	stop = 0
	for k = 1:d
        xk = x + 0
		xk[k] += 1
		fkp = f(xk) - fx
		xk[k] = x[k] - 1
		fkm = f(xk) - fx
		if fkp >= 0 && fkm >= 0
		    stop += 1
        end
		g[k] = (fkp - fkm) * 0.5
	end
	p = g / maximum(abs(g))
	println("g: ", g)
	if stop >= d
        println("no decrease in any direction")
	    return x
	end
	fmin = fx
	xmin = x + 0
    for i = 1:n
	    xk = x - round(Integer, p * i)
		if any( ( xk .< 1 ) | ( xk .> n ) )
		    break
		end
		fxk = f(xk)
			println("line ", i, " val ", fxk, " at ", xk')
		if fxk < fmin
	        fmin = fxk
			xmin = xk
			println("newmin ", fxk, " at ", xk')
		end
	end
	println("minimal value ", fmin, " in ", xmin')
    xmin
end

function opti(x, f, n) 
    while true
	    xa = gradstep(x, f, n)
		if xa == x
		    break
		end
		x = xa + 0
	end
	x
end

function opti!(s::Sums)
  n = length(s.a)
  function f(x)
      deviation(s, [1, x[1], x[2], n])
  end
  x = [n÷3, n÷3*2]
  x = opti(x, f, n)
  s.optind = x
  s.optval = f(x)
  s.optx   = ( s.aa[[x; n]] - s.aa[[1;x]] ) ./ ( [x; n] - [1; x] )
end

