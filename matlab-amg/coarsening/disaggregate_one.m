function [ A_dis, P ] = disaggregate_one( A, i, w )
% disaggreate one vertex
% A is adjacency matrix with 0 diagonal

% size of the orginal graph
n = size(A,1);
is_dis = false(n,1);
is_dis(i) = true;

% compute the disgree of i-th vertice
%d = spones(A)*ones(n,1);
di = nnz(A(i,:));

% size of the new graph
n_dis = (n-1) + di;

% allocate the new graph
A_dis = sparse(n_dis+1,n_dis+1);

% form new graph
% part of graph that we do not touch
%A_dis(1:n-1,1:n-1) = A(~is_dis, ~is_dis);
A_dis(1:n,1:n) = A;
A_dis(i,:) = [];
A_dis(:,i) = [];

% cycle graph for the disaggregated vertex
A_cycle = spdiags([w*ones(di,1), w*ones(di,1)], [-1, 1], di, di);
A_cycle(1,end) = w;
A_cycle(end,1) = w;

A_dis(n:end,n:end) = A_cycle;

% edges between untouched vertices and disaggregated vertex i
[~, col_idx, val] = find(A(i,:));
tmp_idx = (col_idx > i);
col_idx(tmp_idx) = col_idx(tmp_idx) - 1;

A_between = sparse((1:di)', col_idx', val', di, n-1);

A_dis(n:end,1:n-1) = A_between;
A_dis(1:n-1,n:end) = A_between';

% form P
P = sparse(n_dis, n);

P(1:n-1,~is_dis) = speye(n-1,n-1);
P(n:end,i) = ones(di,1);

end

