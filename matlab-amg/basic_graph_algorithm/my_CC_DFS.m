function [ C, num_C, E] = my_CC_DFS( G )
% find all connected components using DFS
%
% @Xiaozhe Hu, Department of Mathematics, Tufts University
% 02/01/2016
%
% Input:       
%       G:      undirected graph
%
% Output: 
%       C:      connected components
%   num_C:      number of connected components
%       E:      egeds on each the BFS tree

%---------------------
% initilize
%---------------------
C = cell(1);         % connected components
E = cell(1);         % DFS edges

k = 0;               % no connected component at beginning          

n = numnodes(G);     % number of vertices

label = -1*ones(n,1);    % label of nodes 
                         % -1: unvisited
                         % 0:  visited
                         % >0: true label
                         
while ~isempty(find(label == -1, 1))
    
    k = k+1;     % we are in the while loop, which means we will
                 % find a new connected component
                 % increase the number of connected component
    
    idx = find(label == -1);    % find unvisited vertex
    [ C{k}, label, E{k} ] = my_MDFS( G, idx(1), label); % run DFS
       
end

num_C = k;

end

