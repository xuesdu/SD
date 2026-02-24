function [x] = forward_gs(A, b, x, DL, nsmooth)
% Forwrd Gauss-Seidel smoother
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % GS iteration
    x = x + DL\(b - A*x);   

end
    
end