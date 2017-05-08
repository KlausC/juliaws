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
    a2 = - a / 2
    if a2 == 0
        dd = -b
		sd = sqrt(abs(b))
    else
        if abs(a2) > abs(b/a2)
            dd = 1 - b / a2 / a2
            sd = abs(a2) * sqrt(1 - b / a2 / a2)
        else
            dd = a2 * a2 - b
            sd = sqrt(abs(dd))
        end
    end
    if dd > 0
        if a < 0
            r1 = a2 + sd
            r2 = b / r1
        else
            r2 = a2 - sd
            r1 = b / r2
        end
    elseif dd == 0
        r1 = a2
        r2 = a2
    else # dd < 0
        r1 = a2 + sd * im
        r2 = conj(r1)
    end
    [r1, r2]
end

# Calculate all roots of a real cubic polynom
# x^3 + a x^2 + b x + c
#
function cubicroots1(a::Real, b::Real, c::Real)
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
	elseif d == 0
		qp3 = q / p * 3 / 2
		z = [ qp3 * 2, -qp3, -qp3 ]
	elseif q == 0
        r1, r2 = squareroots(0, p)
		return [ 0, r1, r2 ]
	else # d != 0
		E = one(d)
		p2 = sqrt(abs(p))
		ff = p2
		p32 = p2 * p2 * p2 
		if p32 <= abs(q)
			ff = cbrt(q)
			pp = p / ff / ff
			qq = E
			x0 = cubic_casea(pp, qq) * ff
		elseif p < 0
			pp = -E
			qq = q / p32
			x0 = cubic_caseb1(pp, qq) * ff
		else
			pp = E
			qq = q / p32
			x0 = cubic_caseb2(pp, qq) * ff
		end
		z = halley(p, q, x0)
	end
	z - a3
end

function cubic_casea(p,q)
	# abs(p) <= 1 and q == 1
	# A = 0.6823278038280193 # A^3 + A - 1 = 0
	# B = 1.324717957244746  # B^3 - A - 1 = 0
    # bb = (A + B - 2) / (B - A)
    # aa = B - 1 - B * bb
	# co = 0.012
	#(aa * p - 1) / (bb * p + 1) - (p*p-1) * p * co

	# rational approximation with x * (a1 + a2*x + a3 *x^2) / (1 + b2*x + b3*x^2)
	# improves initial values extremely -- maximal error 5.3058e-5 on [-1,1]
	const a = [ 0.9999654609142731, -0.6015406943495143, 0.1242305837221518]
	const b = [ 1.0, -0.26842160719590985, 0.03435069390818641]
	num = (a[3] * p + a[2]) * p + a[1]
	nom = (b[3] * p + b[2]) * p + b[1]
	- num / nom
end

function cubic_caseb1(p,q)
	# p == -1 and abs(q) <= 1
	# B = 1.324717957244746  # B^3 - A - 1 = 0
	# co = 0.09
	# co2 = 0.08
    # bb = 3/ 2 - B
    # (bb * q * q - 0.5) * q - sign(q) - (q*q-1) * q * ( co + (q*q-0.22) * co2)

	# rational approximation with x * (a1 + a2*x^2 + a3 *x^4) / (1 + b2*x^2 + b3*x^4)
	# improves initial values extremely -- maximal error 7.43892e-4 on [-1,1]
	const a = [ 0.4696068641627135, 1.5368525961245831, -0.8767730355444439]
	const b = [ 1.0, 4.9999999999975095, -2.5130343845600365]
	q2 = q * q
	num = (a[3] * q2 + a[2]) * q2 + a[1]
	nom = (b[3] * q2 + b[2]) * q2 + b[1]
	- q * num / nom - sign(q)
end

function cubic_caseb2(p,q)
	# p == 1 and abs(q) <= 1
	# A = 0.6823278038280193 # A^3 + A - 1 = 0
	# co = 0.1
	# bb = 1 - A
	# (bb * q * q - 1) * q - (q*q-1) * q * co
	
	# rational approximation with x * (a1 + a2*x^2 + a3 *x^4) / (1 + b2*x^2 + b3*x^4)
	# improves initial values extremely -- maximal error 1.54528e-5i on [-1,1]
	const a = [ 0.9997184848333441, 2.996635923118964, 0.7254540667242799]
	const b = [ 1.0, 3.9820217278022136, 1.9379682167619132]
	q2 = q * q
	num = (a[3] * q2 + a[2]) * q2 + a[1]
	nom = (b[3] * q2 + b[2]) * q2 + b[1]
	- q * num / nom
end


function halley(p, q, x)
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
			# println("x ", x, " ", z1, " ", (x - z1) / z1)
		end
		b0 = z1 * z1 + p 
        r1, r2 = squareroots(z1, b0)
		z = [ z1, r1, r2 ]
	end

function cubicroots(a::Real, b::Real, c::Real)
	E = one(c)
	if abs(c) <= E
		cubicroots1(a,b,c)
    else
        E ./ cubicroots1(b/c, a/c, one(c)/c)
	end
end
trans(t, p) = 1.0 ./ (cbrt(1 + 2 *(p/3).^3) * complex(t) - p/3 )
ptrans(p) = - p^2 / 3 / cbrt(1 + 2 * (p/3)^3)^2

# trans( cubicroots(0, ptrans(p), 1), p ) = cubicroots(0, p, 1)
