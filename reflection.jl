
global types_tree = Dict{DataType, DataType}()

function subtypestree(t::DataType)

  function subtypestree2(t::DataType, ind, parent)
      loop = haskey(D, t) ? "*" : ""
      if doprint
          println(string(ind, t, loop))
      end
      if loop == ""
      D[t] = parent
      s = subtypes(t)
      for t2 in s
          subtypestree2(t2, string(ind, "  "), t)
      end
    end
  end
  
    D = t == Any ? types_tree : Dict{DataType,DataType}()
    doprint = ( t != Any )
    subtypestree2(t, "", Any)
end


function parenttypes(t::DataType)
    if ! haskey(types_tree, t)
	    empty!(types_tree)
        subtypestree(Any)
    end
    types_tree[t]
end



