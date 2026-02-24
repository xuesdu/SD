function [ r ] = SGS_prec( r, DL, DU, D_inv)
% Symmetric Gauss-Seidel preconditioner
%
% @ Junyuan Lin & Xiaozhe Hu, Tufts University


%----------------------------
% Step 1: Main loop
%----------------------------
%for i = 1:nsmooth
    
    % ILU iteration
    r = DU\(D_inv\(DL\r));  

%end

end

