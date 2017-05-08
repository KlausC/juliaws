function fejer2(n::Int)
# Weights of the Clenshaw-Curtis, Fejer1, and Fejer2 quadratures by DFTs
# n>1. Nodes: x_k = cos(k*pi/n)
n2 = mod(n,2);
u0 = 1 / (n^2-1+n2);
# Boundary weights of CC
# Clenshaw-Curtis nodes: k=0,1,...,n; vector of weights wcc, wcc_0=wcc_n=u0
L = 0:n-1;
m = min(L, n-L);
r = 1 ./ (0.5 - 2 * m.^2); # auxiliary vectors for all rules
vc = r - u0;
wcc = real(ifft(vc)); # Clenshaw-Curtis weights
#
# Chebyshev or Fejer-1 nodes: k=1/2,3/2,...,n-1/2; vector of weights wf1
s1 = sign(n/2 - L);
# auxiliary vector
v1 = s1 .* r .* exp(im*pi/n*L);
wf1 = real(ifft(v1)); # Chebyshev or Fejer-1 weights
#
# Filippi or Fejer-2 nodes: k=0,1,...,n; vector of weights wf2, wf2_0=wf2_n=0

flag = abs(n/2 - L) .< 1;
s2 = 1 + flag * (n-(n/2+1)*n2); # auxiliary vectors
v2 = s2 .* r;
wf2 = real(ifft(v2)); # Filippi or Fejer-2 weights
wf2[1] = 0

xf1 = cos(pi/n*((0:n-1)+0.5))
xcc = cos(pi/n*(0:n))
xf2 = xcc
[xf1 wf1], [xf2 [wf2; 0]], [xcc [wcc; wcc[1]]]
end
