function [P] = generate_classical_P_Jacobi(A, isC, radius)
% generate prolongation using direct interpolation (classical AMG)
%
% @ Xiaozhe Hu, Tufts University 

%------------------
% local variables
%------------------
N = size(A,1);
Nc = sum(isC);
Nf = N - Nc;

C_point = find(isC);
F_point = find(~isC);

Aff = A(F_point,F_point);
Afc = A(F_point,C_point);

%------------------
% generate P
%------------------
% initialzie P (fill in the C-point of P)
P = sparse(C_point, (1:Nc)', ones(Nc,1), N, Nc);

% get interpolation for F-points
% initilize W (using aggregation-based AMG)
[~,idx]=max(abs(Afc),[],2);
W = sparse((1:Nf)', idx, ones(Nf,1), Nf, Nc);

% one step of Jacobi to 
Dffinv = spdiags(1./diag(Aff),0,Nf,Nf);
for i=1:radius
    W = W - Dffinv*(Afc+Aff*W);
end

% put W in P
P(F_point,:) = W;

end

