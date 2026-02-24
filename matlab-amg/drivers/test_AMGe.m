
% test AMGe
% 
% @ Xiaozhe Hu, Tufts University
%
% NOTE: to be added to the repository

%---------------------
% generate matrices
%---------------------
tic;
A = assembleLaplace(20);
toc;

%---------------------
% get size
%---------------------
N = size(A,2);

%---------------------
% right hand side
%---------------------
%b = zeros(N,1);
exact_x = rand(N,1) + (1:N)';
b = A*exact_x;

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = 100; % size of the coarest level

amgParam.strong_connection = 0.08; 
amgParam.agg_type = 'HEC';  % 'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening
                                                    %| 'MWM': maximal weighted matching
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve pahse parameters
amgParam.cycle_type = 'V'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle | 'K': K-cycle (nonlinear AMLI-cycle)

amgParam.coarse_it = 2; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K'
amgParam.coarse_krylov_type = 'GCG'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGRMES': flexible GMRes 
                                     
amgParam.smooth_type = 'GS';   % 'Rd': Richardson | 'Jcb': Jacobi | 'GS': Gauss-Seidel | 'K1': Kaczmarz_AAt | 'K12': Kaczmarz_AtA
amgParam.n_presmooth = 1; % number of presmoothing      
amgParam.n_postsmooth =1; % number of postsmoothing
amgParam.smooth_omega = 1.2; % weight for Richardson smoother;

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 100;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-8;    % when AMG is used as standalone solver, tolerance for the reletive residual

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother
% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'nofill';  % nofill, crout, ilutp
amgParam.droptol = 0.1; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = 'CG'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver | 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 1; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 200; % maximal number of iterations that is allowed
iterParam.tol = 1e-6;  % Tolerance for relative residual 
iterParam.restart = 200; % restart number for flexible GMRes

%---------------------
% setup phase
%---------------------
amgData = [];
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
    amgData = AMGe_Setup(A, amgParam);
end

%---------------------
% solve pahse
%---------------------
x = zeros(N,1); %initial guess 

switch iterParam.solver_type
    
    case 'SL', % use simple linear iterative solver
        x = Simple_Solve(A, b, x, iterParam);

    case 'AMG', % use AMG directly
        x = AMG_Solve(amgData, b, x, amgParam);
       
    otherwise, % use preconditioned Krylov methods
        x = Krylov_Solve(A, b, x, iterParam, amgParam, amgData);
        
end

%----------------
% test Af
%----------------
% Af = amgData(1).Af;
% xf = rand(N,1);
% bf = Af*xf;
% 
% iterParam.prec_type = 'NULL';
% x = zeros(N,1);
% 
% xff = Krylov_Solve(Af, bf, x, iterParam, amgParam, amgData);

%-------------------------------
% test strong approximation property
%-------------------------------
% % get Po
% Po = amgData(1).Po;
% 
% % construct coarse grid problem
% Ao = Po'*A*Po;
% fo = Po'*b;
% 
% % get size 
% No = size(Ao, 1);
% 
% % solve Ao*xo = fo
% iterParam.prec_type = 'NULL';
% iterParam.tol = 1e-10;
% 
% %xo = Krylov_Solve(Ao, fo, zeros(No,1), iterParam, amgParam, amgData);
% xo = Ao\fo;
% 
% % get xc = Po*xo;
% xc = Po*xo;
% 
% % check SAP constant
% eta_s = ( (exact_x - xc)'*A*(exact_x - xc) )/( 8*b'*b ) 


