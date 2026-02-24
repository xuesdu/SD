function [ L ] = assembleGraphLaplace(N)
%
% Copyright (C)  Xiaozhe Hu.

e = ones(N,1);
NN = N^2;

L1d = spdiags([-1*e 2*e -1*e], -1:1, N, N);
I = speye(N,N);

L =  kron(L1d, I) + kron(I, L1d);

L = L - spdiags(diag(L), 0, NN, NN);

L = L + spdiags(-sum(L,2), 0, NN, NN);

end

