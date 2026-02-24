function [ z ] = MSP_prec( r, M, H, P_mat, n_step)
% Multilevel Subgraph Preconditioner
%
% @ Junyuan Lin, Loyola Marymount University & Xiaozhe Hu, Tufts University
%n = length(r);
N = size(r,1);
z = zeros(N,1); 
%while norm(r) > tol
tol = 1e-2;
%-------------------
% Main loop 
%-------------------
r_msp = P_mat *r;
z_msp = zeros(length(r_msp),1);
if isempty(H)
    z_msp = M\r_msp;
    
elseif isempty(M)
    z_msp = H\r_msp;
else
    [z_msp, iter, r_msp] = Prec_CG(M, r_msp, z_msp, H, n_step, tol, 0);
end

z = P_mat'*z_msp;
z = z - z' * ones(N, 1) / N;
    
    
    
    %x = pinv(full(H))\r;
    % make r orthogonal to nullspace of A
%     nullspace = null(A);
%     for i = 1:size(nullspace, 2)
%         r = r - r' * nullspace(:,i) / norm(nullspace(:,i))*nullspace(:,i);
%     end
%end

end

