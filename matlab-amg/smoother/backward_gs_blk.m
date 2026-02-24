function [x] = backward_gs_blk(A, b, x, DU_blk, nsmooth)
% Backward block Gauss-Seidel smoother
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % GS iteration
    x = x + DU_blk\(b - A*x);   
     
end