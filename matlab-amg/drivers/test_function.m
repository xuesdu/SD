function test_function(L, H, M, b, x, tol, solver_type, prec_type, threshold, layer_ind, level, P_mat)
% driver for testing Ld
% Inputs: L - graph Laplacian from loading the dataset
%         H - graph Laplacian of tree-structured sparsifier/subgraph
%         M - graph Laplacian of positive subgraph
%         b - right hand side
%         tol - tolerance
%         solver_type, eg: 'AMG'
%         prec_type, eg: 'NULL'
%         threshold - upper bound of avg. convergence rate before restart
%                     input range (0, 1]. 0 means resetup
%                     everytime, and 1 means never change
%                     hierariachies (non-adaptive)
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

%---------------------
% load matrices
%---------------------
% load(filename);
% A = Problem.A;

% % get the largest connected component
% G = graph(A);
% S = conncomp(G);
% i = mode(S);
% largest_conn = find(S==i);
% A = A(largest_conn,:);
% A = A(:,largest_conn);
% 
% % % get a undirected (symmetric) positive weighted graph
% A_sym = (A + A')/2;
% A = abs(A_sym);
% % 
% %take out diagonals of A (if any)
%n = size(H,1);
%H = H;
% d = diag(A);
% A = (A-spdiags(d,0,n,n));
% d = sum(A);
% if size(d,1)==1
%    d = d';
% end
% % 
% % 
% L = spdiags(d,0,n,n)-A;
% % % 
% p = symrcm(L);
% L = L(p,p);

% get normalized graph laplacian
% D_half = spdiags(1./sqrt(d),0,n,n);
% L = D_half * L * D_half;

%N = size(H,1);

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = 1; % size of the coarest level

amgParam.strong_connection = 0.01; 
amgParam.agg_type = 'MWM';  % 'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening 
                                      % | 'path_matching': path cover coarsening
                                      % | 'path_matching_aff': path cover coarsening using affinity
                                      % | 'path_matching_adapt': adaptive path cover coarsening using affinity
                                      % | 'regular_matching': regular matching coarsening
                                      % 'MWM': Maximal Weighted Matching
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve phase parameters
amgParam.cycle_type = 'nV'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle | 'K': K-cycle (nonlinear AMLI-cycle)

amgParam.coarse_it = 3; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K'
amgParam.coarse_krylov_type = 'GCG'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGMRES': flexible GMRes 
                                     
amgParam.smooth_type = 'GS';   % 'Rd': Richardson | 'Jcb': Jacobi | 'GS': Gauss-Seidel | 'K1': Kaczmarz_AAt | 'K12': Kaczmarz_AtA
amgParam.n_presmooth = 1; % number of presmoothing      
amgParam.n_postsmooth =1; % number of postsmoothing
amgParam.smooth_omega = 1.2; % weight for Richardson smoother;

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 2000;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = tol;    % when AMG is used as standalone solver, tolerance for the reletive residual
amgParam.threshold = threshold;    % resetup criterion. Upper bound of avg. convergence rate

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother

% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'nofill';  % nofill, crout, ilutp
amgParam.droptol = 0.01; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.ILU_milu = 'off';  % row, col, off
amgParam.ILU_udiag = 1;
amgParam.ILU_thresh = 1;

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = solver_type; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver 
                              % | 'CG': conjugate gradiant | 'FGMRES':
                              % flexible GMRes | AMG_path: path cover AMG
                              % solver | AMG_Adapt: adaptive AMG solver
                              % 'AMG_path': path cover solver
iterParam.prec_type = prec_type; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 1; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 2000; % maximal number of iterations that is allowed
iterParam.tol = tol;  % Tolerance for relative residual 
iterParam.restart = 100; % restart number for flexible GMRes
% for MSP
iterParam.layer_ind = layer_ind;
iterParam.max_level = level;


% initial guess
%rng(5)
%x = rand(N,1);

%x = zeros(N,1);


%---------------------
% setup phase
%---------------------
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
%if  strcmp(iterParam.prec_type, 'AMG') 
    amgData = AMG_Setup(H, amgParam);
end


%---------------------
% solve phase
%---------------------

%H = L_G; % for solving H^dagger L_G x = H^dagger b

switch iterParam.solver_type
    
    case 'SL', % use simple linear iterative solver
        x = Simple_Solve(L, H, b, x, iterParam);
        
    case 'blk', % use block Gaussian Elimination solver
        x = blkGE_Solve(L, H, b, x, iterParam);

    case 'AMG', % use AMG directly
        x = AMG_Solve(amgData, b, x, amgParam);
    
    case 'AMG_Adapt' % use adaptive AMG solver
        % use exact_x
        x = AMG_Solve_Adapt(L, b, x, amgParam);
        
    otherwise, % use preconditioned Krylov methods
        try
            x = Krylov_Solve_MSP(L, H, M, b, x, iterParam, amgParam, amgData, P_mat);
        catch
            x = Krylov_Solve_MSP(L, H, M, b, x, iterParam, amgParam, [], P_mat);  % when AMG is not called
        end
        %x = x - x' * ones(N,1) / N; % turn this on when solving graph
        %problems
        
end



end

