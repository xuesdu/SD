function [ L ] = assembleGraphLaplace3D(N)
%
% Copyright (C)  Xiaozhe Hu.

e = ones(N,1);
NNN = N^3;

L1d = spdiags([-1*e 2*e -1*e], -1:1, N, N);
I = speye(N,N);

L =  kron(kron(L1d, I),I) + kron(kron(I, L1d),I) + kron(kron(I,I),L1d);

%L = L - spdiags(diag(L), 0, NNN, NNN);

%L = L + spdiags(-sum(L,2), 0, NNN, NNN);

end
