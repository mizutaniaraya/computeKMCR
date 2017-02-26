function C=makeNormalizedRand(d,k)

A=rand(d,k)-1/2;
S=dot(A,A);
C=A*spdiags(sqrt(S').^-1,0,k,k);

return
