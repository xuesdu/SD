function [total_stretch] = calculate_stretch(L_G,L_T)
% A function to calculate the total stretch of spanning tree T of G
% Inputs: L_G - graph Laplacian of G
%             L_T - graph Laplacian of T
% Output: total_stretch (following definition in "A Note on Preconditioning by 
%                                     Low-Stretch Spanning Trees" by
%                                     Spielman and Woo 2009)

% Created by Junyuan Joanne Lin, Loyola Marymount University, Dept. of Math

% initialization
total_stretch = 0;

% find the edges in G but not in T
edges_G = triu(L_G,1);
edges_T = triu(L_T,1);
[row,col] = find(edges_G - edges_T);

%find stretch of each edge
G_T = graph(L_T,'omitself'); % assuming undirected
G_T.Edges.Weight = abs(1./G_T.Edges.Weight);
for i = 1:length(row)
    stretch = 0;
    p = shortestpath(G_T,row(i),col(i)); 
    for j = 1:length(p)-1
        stretch = stretch + 1/abs(L_T(p(j), p(j+1)));
    end
    stretch = stretch*abs(L_G(row(i),col(i)));
    total_stretch = total_stretch + stretch;
end
end

