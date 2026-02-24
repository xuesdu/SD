function [ L_dis, P,  Ds ] = disaggregate_laplacian( L, w, theta )
% form disaggregation
%
% BRIEF: 
%   Disaggregate a given graph (based on graph Laplacian matrix)
%
% INPUT: 
%   L:               Laplacian matrix of original graph
%   w:              weight of disaggregated subgraphs
%   theta:        Percentage of vertices (degree>=3) that will be disaggregated 
%
% OUTPUT:
%   L_dis:        Laplaican matrix of disaggregated graph
%   P:              Prolongtation from original to disaggregated graph
%   Ds:            Diagonal scaling (1/sqrt(degree))   
%
% USAGE:
%    [ L_dis, P,  Ds ] = disaggregate( L, w, theta ): theta<1, some portion of the vertices will be disaggregated
%    [ L_dis, P,  Ds ] = disaggregate( L, w, 1): all the vertices with degree great than 3 will be disaggreted
%  
% COPYRIGHT:
%
% @ Xiaozhe Hu, Tufts University
%
%
% NOTE: need to be added to repository

% get size of original graph
n = size(L,1);

% form adjacency matrix
D = spdiags(diag(L),0,n,n);
A = D-L;

% disaggregated
%[ A_dis, P,  Ds ] = disaggregate( A, w, theta );
[ A_dis, P, Ds ] = disaggregate_adj( A, w, theta );

% form disaggregted graph Laplacian
n_dis = size(A_dis,1);
d_dis = A_dis*ones(n_dis,1);
D_dis = spdiags(d_dis,0,n_dis,n_dis);
L_dis = D_dis - A_dis;

end