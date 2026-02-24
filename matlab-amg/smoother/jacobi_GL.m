function [x] = jacobi_GL(A, D, x, nsmooth)
% Jacobi Smoother for Graph Laplacian
%
% X.Hu & J. Urschel

%----------------------------
% Step 0: Prepare
%----------------------------
n = size(x,1);
m = size(x,2);

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % Jacobi iteration
    for j = 1:m
        x(:,j) = x(:,j) - (A*x(:,j))./D;  
    end
    
    % orthogonalize to the constant vector
    for j = 1:m
        x(:,j) = x(:,j) - sum(x(:,j))/n;
    end
      
    % find orthonormal basis
    [x, ~] = qr(x, 0);
    
end