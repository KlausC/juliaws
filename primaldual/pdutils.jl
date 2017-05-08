#
# julia
#
# Primal-Dual Interior Point Methods for linear programming
# different functions and derivatives to be used in PD-linear algebra
# g(-x) * g(x) = 1 and g(-inf) = 0; continuous derivative.
#
function glin(x::Real)
	e = typeof(x)(1)
    x >= 0 ? (x + e, e) : (e / (e - x), e / (e - x)^2 )
end
@vectorize_1arg Real glin

function gexp(x::Real)
	ex = exp(x)
	(ex, ex)
end
@vectorize_1arg Real gexp

function gsqrt(x::Real)
	e = typeof(x)(1)
	sz = sqrt(x * x + e)
	if x >= 0
		(x + sz, x / sz + e)
	else
		szmz = sz - x
		(e / szmz, e / (szmz * sz))
	end
end
@vectorize_1arg Real gsqrt

#
# Standard potential function
# rho must be > n in to obtain a valid potential
function potential(x, s, rho)
	rho * dot(x, s) - sum(log(x .* s))
end

# Naive way of calculation matrix product A * D * A'
# where A is mxn with m <= n
#
function normalproduct(A, D)

	A * diagm(D) * A'
end

#Basics for legacy PD solvers (ipotential reduction, path-following, Mehrotra)
function pddelta(x, s, D2, A, F, rb, rc, rd)
	
	rbb = A * ( D2 .* rc - rd ./ s ) + rb
	dlam = F \ rbb
	ds = rc - A' * dlam
	dx = - D2 .* ds + rd ./ s

	(dx, ds, dlam)
end

#Basics for legacy PD solvers (potential reduction, path-following, Mehrotra)
function pddelta_aff(x, s, D2, lam, A, F, b, c)
    rb = b - A * x
    rc = c - s - A' * lam
	rd = - x .* s
    printrhs("RA ", rb, rc, rd)
	pddelta(x, s, D2, A, F, rb, rc, rd)
end

#Basics for legacy PD solvers (potential reduction, path-following, Mehrotra)
function pddelta_corr(x, s, D2, lam, A, F, rd)
	
	rb = zeros(lam)
	rc = zeros(x)
    printrhs("RC ", rb, rc, rd)
	pddelta(x, s, D2, A, F, rb, rc, rd)
end

# find maximal steplength alpha in [0,1], so x + alpha * deltax >= 0
function alphamax(x, dx, alphas)
	z = zip(x, dx)
	foldl((a,p) -> p[1] >= 0.0 && p[2] < 0.0 && begin y = -p[1]/p[2]; y < a end ? y : a, alphas, z) 
end

# find maximal x[i]/dx[i] not in 1/(1+beta), 1/(1-beta) for all i
# with beta in (0,1)
function alphamax2(x, dx, beta = 0.1)
    mins, maxs = 0.0, 0.0
    emb, epb = 1.0 - beta, 1.0 + beta
	s = sort(map(p -> p[1] / -p[2], filter(p -> p[2] < 0, zip(x, dx))))
	for si in s
	    sm = si * emb
		sp = si * epb
		if sm > maxs
		    mins, maxs = sm, sp
			if sm > 1.0
			    break
            end
		else
		    maxs = sp
		    if sp >= 1.0
		        break
            end
        end
    end
	maxs >= 1 ? min(mins, 1.0) : 1.0
end

# calcualte main diagonal and normal matrix
function normaldiag(x, s, A)
	D2 = x ./ s
	ADA = normalproduct(A, D2)
	F = factorize(ADA)
	(F, D2)
end

function printpoint(t, x, s, lam)
  if length(x) <= 4
    println(t, "x: ", full(x)')
    println(t, "s: ", full(s)')
    println(t, "λ: ", full(lam)')
  end
end

function printdelta(t, rb, rc, rd)
  if length(rc) <= 4
    println(t, "x: ", norm(rb), " ", full(rb)')
    println(t, "s: ", norm(rc), " ", full(rc)')
    println(t, "λ: ", norm(rd), " ", full(rd)')
  else
    println(t, "x: ", norm(rb))
    println(t, "s: ", norm(rc))
    println(t, "λ: ", norm(rd))
  end
end

function printrhs(t, rb, rc, rd)
  if length(rc) <= 4
    println(t, "rb: ", norm(rb), " ", full(rb)')
    println(t, "rc: ", norm(rc), " ", full(rc)')
    println(t, "rd: ", norm(rd), " ", full(rd)')
  else
    println(t, "rb: ", norm(rb))
    println(t, "rc: ", norm(rc))
    println(t, "rd: ", norm(rd))
  end
end

function printresiduals(t, x, s, lam, sigmue, A, b, c)
  printpoint(t, x, s, lam)
  rb = b - A * x
  rc = c - s - A' * lam
  rd = sigmue - x .* s
  if length(rc) <= 4
    println(t, "rb: ", norm(rb), " ", full(rb)')
    println(t, "rc: ", norm(rc), " ", full(rc)')
    println(t, "rd: ", norm(rd), " ", full(rd)')
  else
    println(t, "rb: ", norm(rb))
    println(t, "rc: ", norm(rc))
    println(t, "rd: ", norm(rd))
  end
end

# count migrations from a old sets to new sets
function migration(old, new)
    m = length(old)
	n = length(new)
	counts = zeros(Int64, m+1, n+1)
	for j in 1:n
	  for i in 1:m
	  	cij = length(intersect(old[i], new[j]))
	  	counts[i,j] = cij
		counts[m+1,j] += cij
		counts[i,n+1] += cij
		counts[m+1,n+1] += cij
      end
	end
	counts
end

# separate variables 1. strong x, 2. strong s, 3. indifferent
function printdistribution(x, s, mue, sar)
	n = length(x)
    so1, so2, so3 = sar[1], sar[2], sar[3]
	s1, s2, s3 = IntSet(), IntSet(), IntSet(1:n)
    for i in 1:n
        xi = x[i]; yi = s[i]
        if xi > mue * yi
            union!(s1, i)
			setdiff!(s3, i)
		elseif yi > mue * xi
            union!(s2, i)
			setdiff!(s3, i)
        end
    end
	
	mig = migration([so1, so2, so3], [s1, s2, s3])
	println("Migration Matrix: \n", mig)
	c1, c2, c3 = length(s1), length(s2), length(s3)
	println("Distribution: ", c1, " ", c2, " ", c3)
	sar[1], sar[2], sar[3] = s1, s2, s3
end

# Infeasible algorithm of Mehrotrai (Algorithm MPC, Wright p198)
function pdmehrotra!(x, s, lam, A, b, c; maxiter = 1000, amax = 0.99, muetol = 1e-6)

  n = length(x)
  mue = dot(x, s) / n
  muemin = mue * muetol
  sar = [IntSet(), IntSet(), IntSet(1:n)]

  for k in 1:maxiter
	println()
	println("Iteration ", k)
    printdistribution(x, s, 10.0, sar)
    printresiduals("", x, s, lam, mue, A, b, c)
	if ( mue <= muemin)
	    break
    end

    F, D2 = normaldiag(x, s, A)
	dax, das, dalam = pddelta_aff(x, s, D2, lam, A, F, b, c)
    printdelta("da", dax, das, dalam)

	apa = alphamax(x, dax, 1.0)
	asa = alphamax(s, das, 1.0)
	xpa = x + apa * dax
	spa = s + asa * das
	lama = lam + asa * dalam
	mua = dot(xpa, spa) / n

	println("apa: ", apa)
	println("asa: ", asa)
    printresiduals("a", xpa, spa, lama, mue, A, b, c)
	println("mue->mua ", mue, "->",  mua)

	sig = min((mua / mue) ^ 3, 1.0)
	println("sig: ", sig)

	rd = sig * mue - dax .* das * apa * asa
	dcx, dcs, dclam = pddelta_corr(x, s, D2, lam, A, F, rd)
    printdelta("dc", dcx, dcs, dclam)

	dx = dcx + dax
	ds = dcs + das
	dlam = dclam + dalam
	ap = min(alphamax(x, dx, Inf) * amax, 1.0)
	as = min(alphamax(s, ds, Inf) * amax, 1.0)
	printdelta("d+", dx, ds, dlam)
	println("ap: ", ap)
	println("as: ", as)

	xc = x + dx * ap
	sc = s + ds * as
	lamc = lam + dlam * as
	muenew = dot(sc, xc) / n
    printresiduals("c", xc, sc, lamc, mue*sig, A, b, c)
	println("mue->munew ", mue, "->", muenew)

	if muenew >= mue * 2.0
		warn("no improvement in mue")
		apa = min(apa, amax)
		asa = min(asa, amax)
		x, s, lam = x + apa * dax, s + asa * das, lama + asa * dalam 
	else
		mue = muenew
		x, s, lam = xc, sc, lamc
	end
  end

end

using MathProgBase
using Clp

# Infeasible algorithm of Mehrotrai (Algorithm MPC, Wright p198)
function pdmehrotra1!(x, s, lam, A, b, c; maxiter = 1000, amax = 0.99, muetol = 1e-6)

  n = length(x)
  mue = dot(x, s) / n
  muemin = mue * muetol
  sar = [IntSet(), IntSet(), IntSet(1:n)]
  solver = ClpSolver(PresolveType=1, SolveType=1, PrimalTolerance=1e-17)

  for k in 1:maxiter
	println()
	println("Iteration ", k)
    printdistribution(x, s, 10.0, sar)
    printresiduals("", x, s, lam, mue, A, b, c)
	if ( mue <= muemin)
	    break
    end

	# Factorising the matrix - main time consuming activity
    F, D2 = normaldiag(x, s, A)
    #
    # Calculate 3 directions for different constraint types
	#

	rb = b - A * x
	rc = c - s - A' * lam
	rd = - x .* s

	dbx, dbs, dblam = pddelta(x, s, D2, A, F, rb, 0.0, 0.0)
	dcx, dcs, dclam = pddelta(x, s, D2, A, F, 0.0, rc, 0.0)
	ddx, dds, ddlam = pddelta(x, s, D2, A, F, 0.0, 0.0, rd)

	# find maximal steplengths ab, ac, ad in [0,1]
	Alin = [dbx dcx ddx]
	blin = -x * amax
	clin = -[norm(rb), norm(rc), norm(rd)] # 
	sol = linprog(clin, Alin, '>', blin, .0, 1.0, solver)
	abx, acx, adx = sol.sol[1], sol.sol[2], sol.sol[3]

	Alin = [dbs dcs dds]
	blin = -s * amax
	clin = -[norm(rb), norm(rc), norm(rd)] # 
	sol = linprog(clin, Alin, '>', blin, .0, 1.0, solver)
	abs, acs, ads = sol.sol[1], sol.sol[2], sol.sol[3]

	xpa = x + [dbx dcx ddx] * [abx, acx, adx]
	spa = s + [dbs dcs dds] * [abs, acs, ads]
	lama = lam + [dblam dclam ddlam] * [abs, acs, ads]
    mua = dot(xpa, spa) / n

	println("ab: ", abx, " ", abs)
	println("ac: ", acx, " ", acs)
	println("ad: ", adx, " ", ads)
	printdelta("da", xpa - x, spa - s, lama - lam)
    printresiduals("", xpa, spa, lama, mue, A, b, c)
	println("mue->mua ", mue, "->",  mua)

	sig = min((mua / mue) ^ 3, 1.0)
	println("sig: ", sig)

	#
	# Corrector step calculation
	#
	rd = sig * mue - (xpa - x) .* (spa - s)
	dcx, dcs, dclam = pddelta_corr(x, s, D2, lam, A, F, rd)
    printdelta("dc", dcx, dcs, dclam)

	dx = dcx + xpa - x
	ds = dcs + spa - s
	dlam = dclam + lama - lam
	ap = min(alphamax(x, dx, Inf) * amax, 1.0)
	as = min(alphamax(s, ds, Inf) * amax, 1.0)
	println("apc: ", ap)
	println("asc: ", as)

	xc = x + dx * ap
	sc = s + ds * as
	lamc = lam + dlam * as
	muenew = dot(sc, xc) / n
    printresiduals("c", xc, sc, lamc, mue*sig, A, b, c)
	println("mue->munew ", mue, "->", muenew)

	if muenew >= mue * 100
		warn("no improvement in mue")
		apa = min(apa, amax)
		asa = min(asa, amax)
		x, s, lam = x + apa * dax, s + asa * das, lama + asa * dalam 
	else
		mue = muenew
		x, s, lam = xc, sc, lamc
	end
  end

end

function driver(f, A, b, c, zeta = 1.0, maxiter = 5)
    x = ones(size(A, 2)) * zeta
    s = ones(size(A, 2)) * zeta
    lam = zeros(size(A, 1))
    f(x, s, lam, A, b, c, maxiter = maxiter, muetol = 1e-17)
end

