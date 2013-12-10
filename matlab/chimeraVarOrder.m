function varOrder = chimeraVarOrder( M, N, L )

indexFlip = M > N;
if indexFlip
  t = M; M = N; N = t;
end

varOrder = zeros(1, M*N*2*L);
vI = 0;

for n=1:N
  for l=1:L
    for m=1:M
      vI = vI + 1;
      varOrder(vI) = chimeraI(m,n,1,l);
    end
  end
end

for n=1:N
  for m=1:M
    for l=1:L
      vI = vI + 1;
      varOrder(vI) = chimeraI(m,n,2,l);
    end
  end
end


  function I = chimeraI(m0, n0, k0, l0)
    if indexFlip 
      I = M*2*L*(n0-1) + 2*L*(m0-1) + L*(2-k0) + l0;
    else
      I = N*2*L*(m0-1) + 2*L*(n0-1) + L*(k0-1) + l0;
    end
  end

end


