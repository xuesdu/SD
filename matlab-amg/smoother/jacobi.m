function [x] = jacobi(A, b, x, Dinv, omega, nsmooth)
% Jacobi Method
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % Jacobi iteration
    x = x + (1/omega)*Dinv.*(b-A*x);          % x_{k+1} = x_k + B(b-Ax_k), B = inv(wD)
    
end