function [A] = Laplacian2Adjacency(L)
%
% BRIEF: 
%   Generate adjacency matrix from the graph Laplacian matrix
%
% INPUT: 
%   L:              graph Laplacian
%
% OUTPUT:
%   A:              adjacency matrix
%
% COPYRIGHT:
% 
%   @ Xiaozhe Hu 10/06/2020, Tufts University

% get size
n = size(L,1);

% get the degree matrix;
D = spdiags(diag(L),0,n,n);

% get the adjacency matrix
A = D - L;

end

