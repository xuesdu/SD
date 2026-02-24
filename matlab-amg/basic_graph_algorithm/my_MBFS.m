function [ C, label, E ] = my_MBFS( G, s, label)
% Breadth First Search (use queue)
%
% @Xiaozhe Hu, Department of Mathematics, Tufts University
% 01/29/2016
%
% Input:       
%       s:      starting vertex
%       G:      undirected graph
%   label:      lable of vertices (-1: unvisited; 0: visited; 1:explored) 
%
% Output: 
%       C:      connected component including vetex s
%       E:      egdes on the BFS tree
%   lable:      updated label of vertices
%

%---------------------
% initilize
%---------------------
n = numnodes(G);                % get number of vertices
A = adjacency(G);               % get adjacency matrix      

C = [];                         % initialize the conncected component
E = [];                         % initialize the edges on the BFS tree

Q = [s];                        % initialize the queue and push s in it
label(s) = 0;                   % label the starting vertex as visited

%---------------------
% main loop
%---------------------
while ~isempty(Q)
   
    v = Q(1);                   % pop the first vetex in Q
    Q(1) = [];
    
    C = [C, v];                 % put v in the connected component
    label(v) = 1;               % mark it as explored
    
    [~,j,~] = find(A(v,:));     % find neighbors of vertex v
    
    for w = j                   % loop over neighbors
        
        if label(w) == -1       % if it is unvisited
            
            label(w) = 0;       % mark it as visited
            Q = [Q, w];         % append it to queue
            E = [E, [v;w]];     % put the edges 
            
        end
        
    end
    
end

end

