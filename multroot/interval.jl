
immutable InfSup{T<:Real}
    inf::T
    sup::T
end

# deliver all powers x^i for i = 1..n
function powers{T}(x::InfSup{T}, n::Integer)
	
	res = Array{InfSup{T}}(n)
    powers!(res, x, n)
end

function powers!{T}(res::Array{InfSup{T}}, x::InfSup{T}, n::Integer)
    if n <= 1
        res[n] = x
    else
        bc = count_ones(-one(n)) - leading_zeros(n-1)
        tp = one(n) << (bc - 1)
        powers!(res, x, tp)
        y = res[tp]
        for i = (tp+1)..n
            res[i] = y * res[i-tp]
        end
    end
    res
end

function *{T}(x::InfSup{T}, y::InfSup{T})
    
