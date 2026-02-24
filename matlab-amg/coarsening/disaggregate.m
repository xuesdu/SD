function [ A_dis, P,  Ds ] = disaggregate( A, w, theta )
% form disaggregation
%
% BRIEF: 
%   Disaggregate a given graph (based on adjacency matrix)
%
% INPUT: 
%   A:              Adjacency matrix of original graph
%   w:              weight of disaggregated subgraphs
%   theta:        Percentage of vertices (degree>=3) that will be disaggregated 
%
% OUTPUT:
%   A_dis:       Adjacency matrix of disaggregated graph
%   P:              Prolongtation from original to disaggregated graph
%   Ds:            Diagonal scaling (1/sqrt(degree))   
%
% USAGE:
%    [ A_dis, P,  Ds ] = disaggregate( A, w, theta ): theta<1, some portion of the vertices will be disaggregated
%    [ A_dis, P,  Ds ] = disaggregate( A, w, 1): all the vertices with degree great than 3 will be disaggreted
%  
% COPYRIGHT:
%
% @ Xiaozhe Hu, Tufts University

% get size
n = size(A,1);

% compute the disgree of each vertex
d = spones(A)*ones(n,1);

% decide which vertices to disaggregate
total_degree = sum(d);
[d_sort, idx] = sort(d, 'descend');

if (theta < 1)
    
    dis_idx = idx(cumsum(d_sort) <= theta*total_degree);
    dis_idx = dis_idx(d(dis_idx) > 3);
    
else
    
     dis_idx = idx(d_sort > 3);
    
end

dis_idx = sort(dis_idx);
num_dis = length(dis_idx);

% prepare for the loop
A_dis = A;
P = speye(n,n);
v_idx = dis_idx;

% main loop
for i = 1:num_dis
     
    % get vertex index
    vertex = v_idx(i);
    
    % disaggregate the vertex
    [A_dis, P_dis] = disaggregate_one(A_dis, vertex, w);
    
     % update P
     P = P_dis*P;
     
     % update indices of vertices
     v_idx = v_idx-1;
    
end

% form diaongal scaling matrix  
n_dis = size(A_dis, 1);
ds = ones(n_dis,1);

start = n-num_dis+1;

for i = 1:num_dis

    % get degree of the vertex
    degree_i = d(dis_idx(i));
    
    % assign scaling
    ds(start:start+degree_i-1) = 1/sqrt(degree_i);
    
    % update start index
    start = start+degree_i;
    
end

Ds = spdiags(ds, 0, n_dis, n_dis);

end

