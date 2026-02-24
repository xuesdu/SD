function [x] = richardson(A, b, x, omega, nsmooth)
% Richardson smoother
%
% @ Xiaozhe Hu, Tufts University

%----------------------------
% Step 1: Main loop
%----------------------------
for i = 1:nsmooth
    
    % GS iteration
    x = x + omega*(b - A*x);

end
    
end
