
# Massepunkt bewegt sich auf Kegelmantel.
# Kegel ist ratationssymmetrisch zu z-Achse.
# Gravitation in Richtung negativer z-Achse: g m/s².
# Öffnungswinkel des Kegels: α rad
# Masse: m kg
# Drehimpuls(z-Komponente) zu Masse: l m²/s
#
# Bewegungsgleichungen aus Lagrangefunktion
# Verallgemeinerte Koordinaten: q = [r; phi].
function Kegel(mi, l, g, α)
	lm = l / m
	sina = sin(α)
	cosa = cos(α)
	function lagrange(q, q1, t)
		r, phi = q
		r1, phi1 = q1
		T = m/2 * (r1^2 + (r*phi1*sina)^2 )
		U = m * g * r * cosa
		T - U
	end
	function cartes(q)
		r, phi = q
		x = r * cos(phi) * cosa
		y = r * sin(phi) * sina
		z = r * cosa
		[x; y; z]
	end
	function mo(t, u)
		r, phi = u[1], u[3]
		r1 = u[2]
		phi1 = lm / (r * sina)^2
		dr = r1
		dphi = phi1
		dr1 = r * (phi1 * sina)^2 - g * cosa
		[dr; dr1; dphi]
	end
	lagrange, cartes, mo
end

