%diameter of the graph
%Input: graph Laplacian
%Output: diameter: diameter of the grpah
% Yue Shen @Tufts University Math Dept.
function diameter=diameter(L)
D = diag(diag(L));
A = D - L;
A = sparse(A);
distances = graphallshortestpaths(A);
[diameter, ~] = max(distances(:));
end