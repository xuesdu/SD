function [ z ] = disagg_exact_prec(r, L_dis, D_s, Ps)
% Disaggregation preconditioner (exact solve on the disaggregate graph)
%
% @ Xiaozhe Hu, Tufts University

% get size
n = size(L_dis,1);

% preconditioning
z =    Ps'*D_s*((L_dis+1e-10*speye(n,n))\(D_s*Ps*r));

end

