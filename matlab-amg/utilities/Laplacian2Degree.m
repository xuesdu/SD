function [degree, D] = Laplacian2Degree(L)
%
% BRIEF: 
%   Generate adjacency matrix from the graph Laplacian matrix
%
% INPUT: 
%   L:              graph Laplacian
%
% OUTPUT:
%   degree:         degree vecrtor
%   D:              degree matrix
%
% COPYRIGHT:
% 
%   @ Xiaozhe Hu 10/06/2020, Tufts University

% get size
n = size(L,1);

% get the degree vector and degree matrix;
degree= diag(L);
D = spdiags(degree,0,n,n);


end