module MPS

export LPData, openmps, readmps!, functions, functions_x2 

#
# object containing all data for a linae problem, which can be derived
# from free form MPS files.
#

type LPData
    name::AbstractString
    objsense::Bool	# true min, false max
	m::Int64     	# number of variables
    n::Int64		# number of constraints

    constnames::AbstractArray{AbstractString,1}
    consttypes::AbstractArray{AbstractString,1}
	varnames::AbstractArray{AbstractString,1}
	
    A::AbstractArray{Float64,2} # problem matrix
    b::AbstractArray{Float64,1} # right hand side
	c::AbstractArray{Float64,1}	# objective function
    
    function LPData()
        ct = Array{AbstractString,1}()
        cn = Array{AbstractString,1}()
        vn = Array{AbstractString,1}()
        A = spzeros(0,0)
        b = zeros(Float64, (0))
        c = zeros(Float64, (0))
        new("", true, 0, 0, ct, cn, vn, A, b, c)
    end
end
#
#
# open mps file
#
function openmps(name::AbstractString)

# try several varaints of file name
# throw exception, is none succeeds
   nam = string(name, ".mps")
   try
       ios = open(string(name, ".mps"))
   catch
       try
           ios = open(string(name, ".MPS"))
       catch
       	   ios = open(name)
       end
   end
end # openmps

function readmps!(data::LPData, ios::IOStream)
    n = 0
    m = 0
    objrow = ""
	
	nop = (x::AbstractString, a::AbstractArray) -> nothing
	detail = nop
    stop = false
	rows = Dict{AbstractString, Int64}()	
	vars = Dict{AbstractString, Int64}()	
	AA = Dict{Tuple{Int64, Int64}, Float64}()
    bb = Dict{Int64, Float64}()
    cc = Dict{Int64, Float64}()

	function objsense(line::AbstractString, la::AbstractArray)
        if la[1] == "MAX" || la[1] == "MAXIMUM"
            data.objsense = false
        elseif la[1] == "MIN" || la[1] == "MINIMUM"
            data.objsense = true
        else
            warn(string("OBJSENSE has invalid value:\"", line, "\"."))
        end
    end

	function processrow(line::AbstractString, la::AbstractArray)
        if la[1] == "N"
            if objrow == ""
                objrow = la[2]
            end
        else
            push!(data.consttypes, la[1])
		    push!(data.constnames, la[2])
            m = length(data.consttypes)
            rows[la[2]] = m
        end
    end

	function processcolumn(line::AbstractString, la::AbstractArray)
        var = la[1]
        if haskey(vars, var)
            k = vars[var]
        else
            n = n + 1
            k = n
            vars[var] = k
        end
        for j in 2:2:length(la)-1
            row = la[j]
            val = float(la[j+1])
            if row == objrow
                cc[k] = val
            elseif haskey(rows, row)
                i = rows[row]
                AA[(i,k)] = val
            end
        end
    end

	function processrhs(line::AbstractString, la::AbstractArray)
        rhs = la[1]
        for j in 2:2:length(la)-1
            i = rows[la[j]]
            val = float(la[j+1])
            bb[i] = val
        end
    end

	function processline(line::AbstractString)
		la = split(line)
    	if line[1] == ' '
            detail(line, la)
        elseif la[1] == "NAME"
            data.name = la[2]
        elseif la[1] == "ROWS"
			detail = processrow
        elseif la[1] == "COLUMNS"
			detail = processcolumn
        elseif la[1] == "RHS"
			detail = processrhs
        elseif la[1] == "RANGES"
			detail = nop
        elseif la[1] == "BOUNDS"
			detail = nop
        elseif la[1] == "SOS"
			detail = nop
        elseif la[1] == "OBJSENSE"
            detail = objsense
        elseif la[1] == "ENDATA" || la[1] == "ENDDATA"
            detail = nop
			stop = true
		else
			warn(string("ignored headline:\"", line, "\"."))
        end    
    end

    for x in eachline(ios)
        processline(chomp(x))
        if stop
            break
        end
    end

	data.m = m
    data.n = n
    data.A = spzeros(m, n)
    data.b = zeros(Float64, (m))
    data.c = zeros(Float64, (n))

	# println(AA)
    # println(bb)
    # println(cc)
    i, k = 0, 0
    v = 0.0
	for kv in AA
        i = kv[1][1]
        k = kv[1][2]
        v = kv[2]
        data.A[i,k] = v
    end
    for kv in bb
        i = kv[1]
        v = kv[2]
        data.b[i] = v
    end
    for kv in cc
        k = kv[1]
        v = kv[2]
        data.c[k] = v
    end
    data.varnames=Array{AbstractString}(n)
    for kv in vars
        data.varnames[kv[2]] = kv[1]
    end
end # readmps!

function functions(p::LPData)
    os = p.objsense ? 1.0 : -1.0
    f(x) = vecdot(x, p.c) * os
	h(x) = p.A * x - p.b
    
	Df(x) = (p.c') * os
	Dh(x) = p.A

    D2L(x, lam) = spzeros(size(p.A,2))

	f, h, Df, Dh, D2L
end # functions

function functions_x2(p::LPData)
    os = p.objsense ? 1.0 : -1.0
    f(x) = vecdot(x.^2, p.c) * 0.5 * os

	h(x) = (p.A * (x.^2) - p.b) * 0.5
    
	Df(x) = (p.c') .* x' * os
	Dh(x) = p.A .* x'

    D2L(x, lam) = p.c' * os + lam * p.A
	f, h, Df, Dh, D2L
end # functions

end # module
