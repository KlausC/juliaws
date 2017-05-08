# Lagrange solver
module Lagrange

function setup(physical::Function, variant, args...)

	lagr, cart, move, startcond, timespan = physical(args...)
	
	problem = ODEProblem(move, startcond(variant), timespan())

	solution = solve(problem)

	xyz = cart.(solution.u)







end
