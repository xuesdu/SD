function [x] = jacobi_blk(A, b, x, Dinv_blk, nsmooth)
% Block Jacobi Method
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 0: prepare
%----------------------------
Dinv = blkdiag(Dinv_blk{:});

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % Jacobi iteration
    x = x + Dinv*(b-A*x);          % x_{k+1} = x_k + B(b-Ax_k), B = inv(D)
    
end