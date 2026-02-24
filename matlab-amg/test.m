
% test AMG: example of using AMG
% 
% @ Xiaozhe Hu, Tufts University

%---------------------
% load matlab-AMG package
%---------------------
setpathAMG;

%---------------------
% generate matrices
%---------------------
% grid-like graphs
N = 256^2;
%N = 64^2;

% usual Laplace on a 2D uniform grid
L = assembleLaplace(sqrt(N));

% Laplace with jumps
%jump_type = 1; % 1: isolated | 2: checker board | 3: random isolated island | 4 random checker board
%bd_type = 1; % 1: Dirichlet | 2: Nuemann
%epsilon = 1e-6; % jump coefficents
%L = assembleJumpLaplace(sqrt(N), epsilon, jump_type, bd_type);

%---------------------
% get size
%---------------------
N = size(L,2);

%---------------------
% right hand side
%---------------------
%exact_x = rand(N,1); %+ (1:N)';
%b = L*exact_x;
b = ones(N,1);
b(N) = 1-N;

%---------------------
% AMG parameters
%---------------------
amgParam = init_AMG_param;

amgParam.print_level = 2; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
                                                    
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = 100; % size of the coarest level

amgParam.amg_type = 'UA';  % 'UA': unsmoothed aggregation AMG 
                           % 'SA': smoothed aggregation AMG
                           % 'C':  classical AMG 

amgParam.strong_connection = 0.0; 
amgParam.agg_type = 'MIS';  % 'HEC': heavy edge coarsening (UA, SA)
                            % 'MIS': maximial indepent set coarsening (UA, SA, C)
                            % 'MWM': maximal weighted matching (UA, SA)
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve pahse parameters
amgParam.cycle_type = 'M'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle |
                            % 'K': K-cycle (nonlinear AMLI-cycle)
                            % 'C': AMLI-cycle using Chebyshev polynomails
                            % 'M': AMLI-cycle using MA
                            % 'N': AMLI-cycle using Nesterov acceleration
                            % 'H': AMLI-cycle using heavy ball
                            % 'Cold': AMLI-cycle using Chebyshev semi
                            % iterative method


amgParam.coarse_it = 2; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K' or 'C' or 'N'
                        % <5 if cycle_type = 'C'
amgParam.coarse_krylov_type = 'GCG'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGRMES': flexible GMRes 
                                     
amgParam.smooth_type = 'GS';   % 'Rd': Richardson | 'Jcb': Jacobi | 'GS': Gauss-Seidel | 'K1': Kaczmarz_AAt | 'K12': Kaczmarz_AtA
amgParam.n_presmooth  = 1; % number of presmoothing      
amgParam.n_postsmooth = 1; % number of postsmoothing
amgParam.smooth_omega = 1.2; % weight for Richardson smoother;

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 300;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-6;    % when AMG is used as standalone solver, tolerance for the reletive residual

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother
% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'crout';  % nofill, crout, ilutp
amgParam.droptol = 0.0001; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = 'CG'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver 
                            % 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 2; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 1000; % maximal number of iterations that is allowed
iterParam.tol = 1e-6;  % Tolerance for relative residual 
iterParam.restart = 100; % restart number for flexible GMRes

%---------------------
% setup phase
%---------------------
amgData = [];
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
    amgData = AMG_Setup(L, amgParam);
end

%---------------------
% solve pahse
%---------------------
rng(1);
x = rand(N,1); %initial guess 
x = x/norm(x);

switch iterParam.solver_type
    
    case 'SL' % use simple linear iterative solver
        [ x, k, err ] = Simple_Solve(L, b, x, iterParam);

    case 'AMG' % use AMG directly
        [ x, k, err ] = AMG_Solve(amgData, b, x, amgParam);
       
    otherwise % use preconditioned Krylov methods
        [ x, k, err ]  = Krylov_Solve(L, b, x, iterParam, amgParam, amgData);
        
end
