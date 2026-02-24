function [x] = forward_gs_blk(A, b, x, DL_blk, nsmooth)
% Forwrd block Gauss-Seidel smoother
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % GS iteration
    x = x + DL_blk\(b - A*x);   

end
    
end