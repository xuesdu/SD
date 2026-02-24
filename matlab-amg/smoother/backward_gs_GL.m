function [x] = backward_gs_GL(A, DU, x, nsmooth)
% Backward Gauss-Seidel Smoother for Graph Laplacian
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
    
    % GS iteration
    for j=1:m
        x(:,j) = x(:,j) - DU\(A*x(:,j));   
    end  
    
    % orthogonalize to the constant vector
    for j=1:m
        x(:,j) = x(:,j) - sum(x(:,j))/n;
    end
        
    % find orthonormal basis
    [x, ~] = qr(x, 0);
    
end