function iterParam = init_Iter_param
% Initialize linear iterative solver parameters
%
% @ Xiaozhe Hu, Tufts University 

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = 'CG'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver | 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 2; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 1000; % maximal number of iterations that is allowed
iterParam.tol = 1e-6;  % Tolerance for relative residual 
iterParam.restart = 1000; % restart number for flexible GMRes

end

