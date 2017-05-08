function fejer(n::Int)
# function [wf1,wf2,wcc] = fejer(n)
# Weights of the Fejer2, Clenshaw-Curtis and Fejer1 quadratures by DFTs
# n>1. Nodes: x_k = cos(k*pi/n)

N=1:2:n-1
l=length(N)
m=n-l
K=0:m-1

# Fejer2 nodes: k=0,1,...,n; weights: wf2, wf2_n=wf2_0=0
v0 = zeros(n+1)
v0[1:l] = 2 ./ N ./ (N-2); v0[l+1] = 1/N[end]
#v0 = [2./N./(N-2); 1/N[end]; zeros(m)]
v2 = -v0[1:end-1] - v0[end:-1:2]
wf2=real(ifft(v2))

# Clenshaw-Curtis nodes: k=0,1,...,n; weights: wcc, wcc_n=wcc_0
g0 = -ones(n)
g0[1+l] += n
g0[1+m] += n
g = g0 / (n^2 - 1 + mod(n,2))
wcc=real(ifft(v2 + g))

# Fejer1 nodes: k=1/2,3/2,...,n-1/2; vector of weights: wf1
a = pi / n * im
v0 = zeros(Complex{Float64}, n+1)
v0[1:m] = exp(a * K) ./ (0.5 - 2.0 * K .^ 2)
#v0 = [2*exp(im*pi*K/n)./(1-4*K.^2); zeros(l+1)]
v1 = v0[1:end-1] + conj(v0[end:-1:2])
wf1=real(ifft(v1))

xf1 = cos(pi/n*((0:n-1)+0.5))
xcc = cos(pi/n*(0:n))
xf2 = xcc
[xf1 wf1], [xf2 [wf2; 0]], [xcc [wcc; wcc[1]]]
end
