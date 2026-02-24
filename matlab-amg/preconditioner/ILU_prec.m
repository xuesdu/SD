function [ r ] = ILU_prec( r, IL, IU)
% ILU preconditioner
%
% @ Junyuan Lin & Xiaozhe Hu, Tufts University


%----------------------------
% Step 1: Main loop
%----------------------------
%for i = 1:nsmooth
    
    % ILU iteration
    r = IU\(IL\r);  

%end

end

