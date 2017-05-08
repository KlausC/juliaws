
# solve the special 2-dimensional convex problem
#
#  maximize{ σ + sum(i)(w[i] * ln(1 - σ * dx[i] -τ * tx[i])) / 0 <= σ <= 1 }
#  
# using Newton's method to solve Lagrange Formula
#

function special2max(dx::AbstractVector, tx::AbstractVector, w::AbstractVector)
  
  n = size(dx,1)
  n == size(tx,1) && n == size(w,1) || error("Differing input vector sizes")
  
  xv(sigma, tau) = 1.0 - sigma * dx - tau * tx
  fs(sigma, tau, x)= -sum(log(x) .* w) - sigma
  function Dfs(sigma, tau, x)
    [ sum(dx ./ x .* w) - 1.0; sum(tx ./ x .* w)]
  end
  function D2fs(sigma, tau, x)
	a11 = sum(dx .* dx .* w ./ x .^ 2)
    a12 = sum(dx .* tx .* w ./ x .^ 2)
	a22 = sum(tx .* tx .* w ./ x .^ 2)
	[a11 a12; a12 a22]
  end
  function maxstep(dsigma, dtau, x::AbstractVector)
    sinv = maximum((dsigma * dx + dtau * tx) ./ x)
	1.0 / max(sinv, 0.5)
  end

  sigma, tau = 0.0, 0.0
  x = xv(sigma, tau)
  fx = fs(sigma, tau, x)
  fy = Inf
  while fx < fy && sigma < 1.0
    y = x; fy = fx
	dsigma, dtau = D2fs(sigma, tau, x) \ Dfs(sigma, tau, x)
	step = maxstep(dsigma, dtau, x)
	println("σ=", sigma, "; τ=", tau, "; dσ=", dsigma, "; dτ=", dtau, "; step=", step)
	step = min(step * 0.5, 1.0)
	sigma -= step * dsigma
	tau -= step * dtau
	sigma = clamp(sigma, 0.0, 1.0)
	x = xv(sigma, tau)
	fx = fs(sigma, tau, x)
  end
  sigma, tau 
end

