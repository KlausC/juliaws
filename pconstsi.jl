"""
Physical constants 
"""
module PhysicalConstantsSI

using SIUnits
using SIUnits.ShortUnits
export Wb, T, H

# unit symbols missing in SIUnits.ShortUnits
One = m / m				# dimensionless
Wb = V * s				# weber - magnetuc flux
T = Wb / m^2			# tesla - magnetic flux density
H = Wb / A				# henry - inductance

import Base: show,*,/,^,+,-

immutable PhysicalConstant{N <: Real,m,kg,s,A,K,mol,cd,rad,sr}
	value::N
	abstol::N
	physdim::SIUnits.SIUnit{m,kg,s,A,K,mol,cd,rad,sr}
	name::AbstractString
	description::AbstractString
end

typealias PC2{N2 <: Real,m2,kg2,s2,A2,K2,mol2,cd2,rad2,sr2} PhysicalConstant{N2,m2,kg2,s2,A2,K2,mol2,cd2,rad2,sr2}

macro def(c, v, tol, ty, d)
	c1 = string(c)
	esc(:(export $c; $c = pu($v, $tol, $ty, $c1, $d)))
end

function pu{N <: Real,NT <:Real,m,kg,s,A,K,mol,cd,rad,sr}(val::N, tol::NT, un::SIUnits.SIUnit{m,kg,s,A,K,mol,cd,rad,sr}, name::AbstractString, description::AbstractString)

	PhysicalConstant(val, N(tol), un, name, description)
end

function show(io::IO, c::PhysicalConstant)
	print(io, "$(c.name): $(c.description): $(c.value) ± $(c.abstol) $(c.physdim)")
end

(+)(a::PhysicalConstant, b::PhysicalConstant) = pu(a.value+b.value, a.abstol+b.abstol, a.physdim, "", "")
(-)(a::PhysicalConstant, b::PhysicalConstant) = pu(a.value-b.value, a.abstol+b.abstol, a.physdim, "", "")
(*){M <: Real}(a::PhysicalConstant, b::M) = pu(a.value*b, a.abstol*b, a.physdim, "", "")
(*){M <: Real}(a::M, b::PhysicalConstant) = pu(a*b.value, a*b.abstol, b.physdim, "", "")
(*)(a::PhysicalConstant, b::PC2) = pu(a.value*b.value, a.abstol*abs(b.value)+abs(a.value)*b.abstol, a.physdim * b.physdim, "", "")
(/){M <: Real}(a::PhysicalConstant, b::M) = pu(a.value/b, a.abstol*b, a.physdim, "", "")
(/){M <: Real}(a::M, b::PhysicalConstant) = pu(a/b.value, a*b.abstol/b.value^2, b.physdim^-1, "", "")
(/)(a::PhysicalConstant, b::PC2) = pu(a.value/b.value, (a.abstol*abs(b.value)+abs(a.value)*b.abstol)/b.value^2, a.physdim / b.physdim, "", "")



# Universal
@def c 		299792458.0		 	0						m/s		"speed of light in vacuum"
@def ν_cs133	9192631770.0	0						Hz		"frequency of radiation corresponding to the transition between the two hyperfine levels of the ground state of the caesium-133 atom"
@def h		6.626070040e-34		0.000000081e-34			J * s	"Planck constant"
@def G		6.67408e-11			0.00031e-11			N/kg^2/s^2	"Newtonian constant of gravitation"

# electromagnetic
@def e0		1.6021766208e-19	0.0000000098e-19		C		"elementary charge"

# physico-chemical
@def NA		6.022140857e23		0.000000074e23			mol^-1	"Avogadro constant"
@def mu		1.660539040e-27		0.000000020e-27			kg		"atomic mass constant 1/12 C12"
@def k		1.38064852e-23		0.00000079e-23			J/K		"Boltzmann constant"
@def p_atm	101325.0			0						Pa		"standard atmosphere"

# atomic and nuclear 
@def R∞		10973731.568508		0.000065				m^-1	"Rydberg constant"	

# electron
@def me		9.10938356e-31		0.00000011e-31			kg		"electron mass"
@def μe		-928.4764620e-26	0.0000057e-26			J / T	"electron magnetic moment"
# muon
@def mμ		1.883531594e-28		0.000000048e-31			kg		"muon masss"
@def μμ		-4.49044826e-26		0.00000010e-26			J/T		"muon magnetic moment"
# tau minus
@def mτ		3.16747e-27			0.00029e-27				kg		"tau mass"
# proton
@def mp		1.672621898e-27		0.000000021e-27			kg		"proton mass"
@def μp		1.4106067873e-26	0.0000000097e-26		J/T		"proton magnetic moment"
@def rp		0.8751e-15			0.0061e-15				m		"proton rms charge radius"
# neutron
@def mn		1.674927471e-27		0.000000021e-27			kg		"neutron mass"
@def μn		-0.96623650e-26		0.0000023e-26			J/T		"neutron magnetic moment"
# deuteron
@def md		3.343583719e-27		0.000000041e-27			kg		"deuteron mass"
@def μd		0.4330735040e-26	0.0000000036e-26		J/T		"deuteron magnetic moment"
@def rd		2.1413e-15			0.0025e-15				m		"deuteron rms charge radius"
# triton
@def mt		5.007356665e-27		0.000000062e-27			kg		"triton mass"
@def μt		1.504609503e-26		0.000000012e-26			J/T		"triton magnetic moment"
# helion
@def mh		5.006412700e-27		0.000000062e-27			kg		"helion mass"
@def μh		-1.074617522e-26	0.000000014e-26			J/T		"helion magnetic moment"
# alpha particle
@def mα		6.644657230e-27		0.000000082e-27			kg		"alpha particle mass"

# aliases and derived constants
# universal
@def μ0	4pi*1e-7				0						N/A^2	"magnetic constant" 
@def ϵ0	1/(μ0*c^2).value		0						F/m		"electric constant"
@def Z0	(μ0*c).value			0		μ0.physdim * c.physdim	"characteristic impedance of vacuum"
@def ℏ	h.value/(2pi)			h.abstol/(2pi)			J * s	"Planck constant slash"

# electromagnetic
@def e_h	e0.value/h.value	0.000000015e14			A/J		"e / h"
@def Φ0		h.value/e0.value/2	0.000000013e-15			J/A		"magnetic flux quantum"
@def G0	2*e0.value^2/h.value	0.0000000018e-5			S		"conductance quantum"
@def G0_1	1/G0.value			0.0000029e4				Ohm		"inverse of conductance quantum"
@def KJ		1/Φ0.value			0.0030e9				Hz/V	"Josephson constant"
@def RK		2/G0.value			0.0000059				Ohm		"von Klitzig constant"

# atomic and nuclear 
@def α	e0.value^2/(ϵ0.value*h.value*c.value*2)	1.7e-12	One		"fine-structure constant"
@def α_1	1/α.value			0.000000031				One		"inverse fine-structure constant"
@def μB		(e0*ℏ/2me).value	0						m^2*A	"Bohr magneton"
# physico-chemical
@def R		(k*NA).value		0.0000048				J/mol/K	"molar gas constant"
@def σ	(pi^2*k^4/(60ℏ^3*c^2)).value	0			kg/s^3/K^4	"Stefan-Boltzmann constant"
end # module
