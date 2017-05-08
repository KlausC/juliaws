
function hilbert(n::Integer, m::Integer = n, S::Type = Float64)
    h = zeros(S,n,m)
    for i = 1:n for j = 1:m
        h[i,j] = S(i+j-1)
    end end
    1 ./ h
end

