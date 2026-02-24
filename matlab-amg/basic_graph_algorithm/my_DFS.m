function [ C, E ] = my_DFS( G, s )
% Depth First Search (use stack)
%
% @Xiaozhe Hu, Department of Mathematics, Tufts University
% 01/29/2016
%
% Input:       
%       s:      starting vertex
%       G:      undirected graph
%
% Output: 
%       C:      connected component including vetex s
%       E:      egeds on the DFS tree
%


%---------------------
% initilize
%---------------------
n = numnodes(G);         % get number of vertices

A = adjacency(G);        % get adjacency matrix

C = [s];                 % connected component
E = [];                  % edges on DFS tree

label = false(n,1);      % label of nodes

S = [s];                 % initialize stack
label(s) = true;         % label the starting vertex as visited

flag = true;

while ~isempty(S)
   
    while flag
        
        v = S(1);                           % get the first vertex in stack
        [~,Nv,~] = find(A(v,:));            % find its neighbors
        idx = find(label(Nv) == false);     % check any adjacent vertex is unvisited
        
        if ~isempty(idx)                    % if there is such vertex
            
            w = Nv(idx(1));                 % choose one of such vertex
            label(w) = true;                 % label it us visited
            C = [C,w];                      % put it in connected component
            E = [E, [v;w]];                 % put the edges 
            S = [w, S];                     % push onto stack
        
        else
            
            flag = false;                   % otherwise, get out of the loop
            
        end
        
    end
    
    S(1) = [];                              % pop stack (backtrack)
    flag = true;                            
        
end