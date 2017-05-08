# This macro defines a global constant variable, but does not pollute the namespace of 
# the current module
# Due to KristofferC - feature request issue 15056 (static local variables missing in julia)

macro static(init)
  var = gensym()
  eval(current_module(), :(const $var = $init))
  var = esc(var)
  quote
    global $var
    $var
  end
end
#
# usage:
#function foo()
#    J = @static zeros(5,5)
#end
#
#@allocated foo()
