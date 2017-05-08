
function normupv(u, v)
  norm(u + v)
end

function normuv(u, v)
  norm(u .* v)
end

function normquot(u, v)
  normuv(u, v) / normupv(u, v) ^2
end

function randquot(n::Int)
  u = rand(n) - 0.5
  v = rand(n) - 0.5
  if vecdot(u, v) < 0
    v = -v
  end
  sc = hypot(norm(u), norm(v))
  u = u / sc
  v = v / sc
  normquot(u, v), u, v
end

