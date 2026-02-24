% generate graph Laplacian on NxN grid
N = 3;
L = assembleGraphLaplace(N);

% get graph 
G = graph(-L, "omitselfloops");

% get incidence matrix
E = incidence(G)';

% get edge weight matrix;
w = G.Edges.Weight;
m = length(w);
W = spdiags(w,0,m,m);


