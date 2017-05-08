
using Base.SparseMatrix.SPQR: solve, qmult, Factorization, QX, QTX, RX_EQUALS_B, RETX_EQUALS_B
using Base.SparseMatrix.CHOLMOD: VTypes, Dense
using Base.SparseMatrix: droptol! 

import Base.getindex

function getindex{T<:VTypes}(qr::Factorization{T}, part::Symbol)
  m, n = size(qr)  
  if part == :Q
    sparse(qmult(QX, qr, Dense(eye(m,n))))
  elseif part == :RI
    sparse(solve(RX_EQUALS_B, qr, Dense(eye(m,n))))
  elseif part in ( :E, :EP )
    ERI = solve(RETX_EQUALS_B, qr, Dense(eye(m,n)))
    EP = sortperm(map( k -> mod(findfirst(ERI[k,:])+n, n+1)+1, 1:n))
    EP = sortperm(EP)  # inverse permutation
    if part == :E
      sparse(1:n, EP, 1.0)
    else
        EP
    end
  end
end
