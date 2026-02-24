function [x, iter] = testAMG(A, b, tol)

N = size(A,1);

% make it SPD
%[~, col] = find(A(1,:));
%A= A + sparse([1;col(2);1;col(2)], [1;1;col(2);col(2)], ones(4,1), N, N);

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 5; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 2; % maximal number of level in AMG
amgParam.coarsest_size = 100; % size of the coarest level

amgParam.amg_type = 'UA';  % 'UA': unsmoothed aggregation AMG 
                           % 'SA': smoothed aggregation AMG
                           % 'C':  classical AMG 

amgParam.strong_connection = 0.0; 
amgParam.agg_type = 'HEC'; %  'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening
                                                    %| 'MWM': maximal weighted matchingamgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve pahse parameters
amgParam.cycle_type = 'V'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle | 'K': K-cycle (nonlinear AMLI-cycle)

amgParam.coarse_it = 2; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K'
amgParam.coarse_krylov_type = 'FGRMES'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGRMES': flexible GMRes 
                                     
amgParam.smooth_type = 'GS';   % 'Rd': Richardson | 'Jb': Jacobi | 'GS': Gauss-Seidel
amgParam.n_presmooth = 1; % number of presmoothing      
amgParam.n_postsmooth = 1; % number of postsmoothing
amgParam.smooth_omega = 1.2; % weight for Richardson smoother;

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 10;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-10;    % when AMG is used as standalone solver, tolerance for the reletive residual

% parameters for Schwarz smoother
amgParam.Schwarz_level = 0;  % how many levels use Schwarz smoother (from fine to coarse), 0 means no Schwarz smoother

% parameters for ILU smoother
amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother

% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'crout';  % nofill, crout, ilutp
amgParam.droptol = 0.0001; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

amgParam.number_eigen = 0;

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = 'FGMRES'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver | 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 5; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 100; % maximal number of iterations that is allowed
iterParam.tol = tol;  % Tolerance for relative residual 
iterParam.restart = 100; % restart number for flexible GMRes

%---------------------
% setup phase
%---------------------
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
    amgData = AMG_Setup(A, amgParam);
end

%---------------------
% solve pahse
%---------------------
x = zeros(N,1); %initial guess 

switch iterParam.solver_type
    
    case 'SL' % use simple linear iterative solver
        [x, iter] = Simple_Solve(A, b, x, iterParam);

    case 'AMG' % use AMG directly
        [x, iter] = AMG_Solve(amgData, b, x, amgParam);

    otherwise % use preconditioned Krylov methods
        [x, iter] = Krylov_Solve(A, b, x, iterParam, amgParam, amgData);
        
end


end

