function [ comp,conn_P,conn_d,conn_A] = get_conncomp( A,P,d )
%check if graph of weighted adjacency matrix is connected. Take the largest
%component for examination (not working for directed graph yet)

% Junyuan Lin and Xiaozhe Hu @Tufts University Math Dept.

%   Input: A-weighted adjacency matrix
%   Output: P-state transition matrix of the largest component
%           d-degree array of the largest component

G = graph(A);
S = conncomp(G);

num = union(S,S); 
num = length(num);%number of connected components
comp = cell(num,1);
conn_P = cell(num,1);
conn_A = cell(num,1);
conn_d = cell(num,1);

for i = 1:num
    comp{i} = find(S==i);
    conn_comp = P(comp{i},:);
    conn_comp = conn_comp(:,comp{i});
    conn_P{i} = conn_comp;
    conn_d{i} = d(comp{i});
    conn_comp = A(comp{i},:);
    conn_comp = conn_comp(:,comp{i});
    conn_A{i} = conn_comp;
end



end

