function [ A_dis, P, Ds ] = disaggregate_adj( A, w, theta)
% form disaggregation
%
% BRIEF: 
%   Disaggregate a given graph (based on adjacency matrix)
%
% INPUT: 
%   A:              Adjacency matrix of original graph
%   w:              weight of disaggregated subgraphs
%   theta:        Percentage of vertices (degree>average degree) that will be disaggregated 
%
% OUTPUT:
%   A_dis:       Adjacency matrix of disaggregated graph
%   P:              Prolongtation from original to disaggregated graph
%   Ds:            Diagonal scaling (1/sqrt(degree))   
%
% USAGE:
%    [ A_dis, P,  Ds ] = disaggregate_adj( A, w, theta): theta<1, some portion of the vertices will be disaggregated
%    [ A_dis, P,  Ds ] = disaggregate_adj( A, w, 1): all the vertices with degree great than average degree will be disaggreted
%  
% COPYRIGHT:
%
% @ Xiaozhe Hu, Tufts University
%
%
% NOTE: need to be added to repository

% get size
n = size(A,1);
[row_A, ~, weight_A] = find(A);

% compute the disgree of each vertex
d = spones(A)*ones(n,1);

% decide which vertices to disaggregate
total_degree = sum(d);
ave_degree = ceil(total_degree/n);
[d_sort, idx] = sort(d, 'descend');

if (theta < 1)
    
    dis_idx = idx(cumsum(d_sort) <= theta*total_degree);
    dis_idx = dis_idx(d(dis_idx) > (max(ave_degree,10)) );
    
else
    
     dis_idx = idx(d_sort > (max(ave_degree, 10)));
    
end

dis_idx = sort(dis_idx);
num_dis = length(dis_idx);
dis_flag = false(n,1);
dis_flag(dis_idx) = true;

% prepare for the disaggregated graph
n_dis = n - num_dis + sum(d(dis_idx));

% aggregation
aggregation = zeros(n_dis, 1);
aggregation(1:n) = (1:n)';

% main loop
local_start = n+1;

vtx_idx = [];

row_C = [];
col_C = [];
weight_C = [];

for vertex = 1:n
     
    if dis_flag(vertex) == true  % current vertex needs to be disaggregated
    
        % local index
        local_end = local_start + d(vertex)-2;
        local_idx = [vertex, local_start:local_end];

        % update aggregation
        aggregation(local_start:local_end) = vertex;

        % update vertex index
        vtx_idx = [vtx_idx; local_idx'];
 
        % add local cycle
        row_C = [row_C; vertex; (local_start:local_end)'; (local_start:local_end)'; vertex];
        col_C = [col_C; (local_start:local_end)'; vertex; vertex; (local_start:local_end)'];
        weight_C = [weight_C; w*ones(d(vertex),1);w*ones(d(vertex),1)];

        % update local start
        local_start = local_end + 1;
    
    else
        
        vtx_idx = [vtx_idx; vertex*ones(d(vertex),1)];
        
    end
    
end

col_A = vtx_idx;
[~, row_idx] = sort(row_A);
row_A(row_idx) = col_A;

A_dis = sparse([row_A; row_C], [col_A; col_C], [weight_A; weight_C], n_dis, n_dis);

P = sparse((1:n_dis)', aggregation,ones(n_dis,1), n_dis, n);

ds = 1./sqrt(d(aggregation));
ds(~dis_flag) = 1; 
Ds = spdiags(ds, 0, n_dis, n_dis);

end

