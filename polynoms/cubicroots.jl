# Calculate all roots of real quadratic polynom
# x^2 + a x + b
#
function squareroots(a::Real, b::Real)
	E = one(b)
	if abs(b) <= E
		squareroots1(a, b)
	else
		E ./ squareroots1(a/b, E/b)
	end
end

function squareroots1(a::Real, b::Real)
	dd = a * a / 4 - b
    sd = sqrt(abs(dd))
    if dd >= 0
        if a < 0
            r1 = -a / 2 + sd
            r2 = b / r1
        else
            r2 = -a / 2 - sd
            r1 = b / r2
        end
    else
        r1 = - a/2 + sd * im
        r2 = conj(r1)
    end
    [r1, r2]
end

# Calculate all roots of a real cubic polynom
# x^3 + a x^2 + b x + c
#
function cubicroots1a(a::Real, b::Real, c::Real)
	a3 = a / 3
	aa3 = a3 * a
    p = b - aa3
    q = ( aa3 * 2 / 3 - b ) * a3 + c
	d = (q / 2) ^ 2 + (p / 3) ^ 3 # discriminant
	if c == 0
        r1, r2 = squareroots(a, b)
		return [ 0, r1, r2 ]
	elseif q == 0 && p == 0
		z = [ 0, 0, 0 ]
	elseif q == 0
        r1, r2 = squareroots(0, p)
		return [ 0, r1, r2 ]
	elseif d >= 0
		sd = sqrt(d)
		if q < 0
			u = cbrt( -q / 2 + sd )
			v = - p / 3 / u
		else
			v = cbrt( -q / 2 - sd )
			u = - p / 3 / v
		end
		upv2 = (u + v) / 2
		umv3 = (u - v) / 2 * sqrt(3)
		z = [ u + v, -upv2 + umv3 * im, -upv2 - umv3 * im ]
	elseif d == 0
		qp3 = q / p * 3 / 2
		z = [ qp3 * 2, -qp3, -qp3 ]
	else # d < 0 && q != 0 - then also p < 0
		x = sqrt(-p) * sign(-q)
		x += q / p * 0.4 
		f(x) = ( x * x + p ) * x + q
	    df(x) = x * x * 3 + p
		d2f(x) = x * 6
		function halley(x)
			fx = f(x)
			dfx = df(x)
			d2fx = d2f(x)
		   # x - fx / dfx # newton
		   # x - fx * dfx / ( dfx * dfx - fx * d2fx / 2 ) # halley
			x - fx / dfx / ( 1 - fx * d2fx / dfx / dfx / 2 ) # halley
           # ( x*dfx*dfx - fx*(x*d2fx/2 +dfx)) / (dfx *dfx - fx * d2fx / 2) # halley
		end
        z1, x, y = x, NaN, NaN
        while z1 != x && z1 != y
			y = x
			x = z1
			z1 = halley(z1)
			println("x ", x, " ", z1, " ", (x - z1) / z1)
		end
		b0 = z1 * z1 + p 
        r1, r2 = squareroots(z1, b0)
		z = [ z1, r1, r2 ]
	end
	z - a3
end

function cubicrootsa(a::Real, b::Real, c::Real)
	E = one(c)
	if abs(c) <= E
		cubicroots1a(a,b,c)
    else
        E ./ cubicroots1a(b/c, a/c, one(c)/c)
	end
end

