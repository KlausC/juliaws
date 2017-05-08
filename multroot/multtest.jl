#
module MultTest

export X
export ptest43, ptest51, ptest52, ptest53


import Polynomials
using Polynomials

# Standard polynomial
X = poly([0.0])
#some test polynomials used by Zeng

# polynomial with roots 1:4 to the power of k
# k = 1..8
function ptest46(k::Int)
    ((X-1)^4 * (X-2)^3 * (X-3)^2 * (X-4) ) ^k
end

# polynomial with non-integer roots
function ptest51(dig::Int)
   p = (X - 10/11)^5 * (X - 20/11)^5 * (X - 30/11)^5
   Poly(map(i->round(p.a[i], dig), 1:length(p.a)))
end

# polynomial with adjacent roots close to 1 with given multiplicity
# epsilon = 0.1, 0.01, 0.001 ...
function ptest52(epsilon::Number)
    (X - 1 + epsilon)^20 * (X - 1)^20 * (X + 0.5)^5
end

# polynomial with given roots (including all conjugate roots) to the power of 2^k
function ptest53(k::Int)
	c = [0.5+im, -1.0+0.2im, -0.1+im, +0.8+0.6im, -0.7+0.7im,
		 1.4, -0.4+0.9im, 0.9, -0.8+0.3im, 0.3+0.8im, 0.6+0.4im]
 
    fX = Poly([1.0])
    for i = 1:length(c)
        a = c[i]
        if isreal(a)
            fX = fX * (X - a)
        else
            fX = fX * ((X-real(a))^2 + imag(a)^2)
        end
    end
	
	for i = 1:k
        fX = fX * fX
    end
    fX
end

end

import MultTest
using MultTest
