
# Doppelpendel in einer Ebene x, z. 
#
# Erster Punkt Masse m1, Länge l1 rotiert um Ursprung.
# Zweiter Punkt Masse m2, Länge l2 rotiert um ersten Punkt.
# Gravitation in Richtung der z-Achse.
#
# Verallgemeinerte Koordinaten: q = [phi1; phi2].
# phi1, phi2 Abweichungen der Pendel zur Senkrechten.
#
function Doppelpendel(m1, l1, m2, l2, g, gamma = 0.0)
	function lagrange(q, q1, t)
		phi1, phi2 = q
		phi1t, phi2t = q1
		T = (m1+m2)/2 * (l1*phi1t)^2 + im2 *((l2*phi2t)^2 + 2l1*l2 *cos(phi1-phi2) * phi1t * phi2t)
		U = -(m1+m2) * g * l1 * cos(phi1) - m2 * g * l2 * cos(phi2)
		T - U
	end
	function raleigh(q, qt)
		phi1t, phi2t = qt
		1/2 * gamma .* [phi1t.^2, (phi1t-phi2t)^2]
	end

	function cartes(q)::Vector
		phi1, phi2  = q
		x1, y1 = l1 * sin(phi1), -l1 * cos(phi1)
		x2, y2 = x1 + l2 * sin(phi2), y1 - l2 * cos(phi2)
		[x1; y1; x2; y2]
	end
	function mo(t, u)::Vector
		phi1, phi2, phi1t, phi2t = u
		dphi1, dphi2 = phi1t, phi2t
		cosdel = cos(phi1-phi2)
		sindel = sin(phi1-phi2)
		m2l1l2cd = m2 * l1 * l2 * cosdel
		M = [(m1+m2)*l1^2 m2l1l2cd; m2l1l2cd m2*l2^2]
		b = [phi1t^2; phi2t^2] * m2 * l1 * l2 * sindel - [(m1+m2)*l1*sin(phi1); m2*l2*sin(phi2)]*g
		b -= gamma .* [ phi1t; phi2t-phi1t] + [gamma[2] * (phi1t-phi2t); zero(phi2t)] 
		println("$t $(det(M))")
		[ dphi1; dphi2; M \ b]
	end
	lagrange, cartes, mo
end

