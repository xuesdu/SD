
% test AMG with elimination
% 
% @ Xiaozhe Hu, Tufts University

%---------------------
% generate matrices
%---------------------
%tic;
% grid-like graphs
% N = 64^2;
% L = assembleGraphLaplace(sqrt(N));

% E-R random graphs
%N = 4000;
%p = 2*log(N)/N;
%[Adj, G] = Erdos_Reyni_Random_Graph( N, p );
%L = laplacian(G); % N^2

% B = incidence(G);
% Bt = B';
% L = spdiags(sum(Adj,2),0, N, N) - Adj;

% W-S random graphs
%N = 2^10;   % nodes
%K = 5;     % average degree (2*K)
%beta = 0;  % beta = 0 is a ring lattice, and beta = 1 is a random graph
%[~, L] = WattsStrogatz(N, K, beta); % N^2


% B-A random graphs
% N = 2000;
% A = ba_net('N', N);
% G =graph(A, 'OmitSelfLoops');
% L = laplacian(G);

%base = 7;
%L = make_poorly_conditioned(L, base);

%toc;

%---------------------
% get size
%---------------------
N = size(L,2);

%---------------------
% right hand side
%---------------------
%b = zeros(N,1);
%exact_x = rand(N,1) + (1:N)';
%b = L*exact_x;
b =zeros(N,1);
b(2) = 1;
b(3) = -1;
%b = buv;

%b = b - sum(b)/N;

%---------------------
% elimination
%---------------------
d = diag(L);
el_idx = (d<=2);
s_idx = (~el_idx);

L11 = L(el_idx, el_idx);
L12 = L(el_idx, s_idx);
L21 = L(s_idx, el_idx);
L22 = L(s_idx, s_idx);

b1 = b(el_idx);
b2 = b(s_idx);

elimination_start = tic;
% form schur complement
[L11L, L11U] = lu(L11);

Ls = L22 - L21*( L11U\(L11L\L12) );

% form right hand side
bs = b2 - L21*(L11\b1);

% get size
Ns = size(Ls,1);

elimination_duration = toc(elimination_start);

fprintf('----------------------------------------------------\n');
fprintf('        elimination costs %f seconds\n', elimination_duration);
fprintf('----------------------------------------------------\n');

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 2; % maximal number of level in AMG
amgParam.coarsest_size = 100; % size of the coarest level

amgParam.strong_connection = 0.00; 
amgParam.agg_type = 'HEC';  % 'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening
                                                    %| 'MWM': maximal weighted matching
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve pahse parameters
amgParam.cycle_type = 'V'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle | 'K': K-cycle (nonlinear AMLI-cycle)

amgParam.coarse_it = 2; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K'
amgParam.coarse_krylov_type = 'GCG'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGRMES': flexible GMRes 
                                     
amgParam.smooth_type = 'Jcb';   % 'Rd': Richardson | 'Jcb': Jacobi | 'GS': Gauss-Seidel | 'K1': Kaczmarz_AAt | 'K12': Kaczmarz_AtA
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
iterParam.print_level = 0; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 3000; % maximal number of iterations that is allowed
iterParam.tol = 1e-8;  % Tolerance for relative residual 
iterParam.restart = 1000; % restart number for flexible GMRes

%---------------------
% setup phase
%---------------------
amgData = [];
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
    amgData = AMG_Setup(Ls, amgParam);
end

%---------------------
% solve pahse
%---------------------
xs = zeros(Ns,1); %initial guess 

switch iterParam.solver_type
    
    case 'SL', % use simple linear iterative solver
        xs = Simple_Solve(Ls, bs, xs, iterParam);

    case 'AMG', % use AMG directly
        xs = AMG_Solve(amgData, bs, xs, amgParam);

    otherwise, % use preconditioned Krylov methods
        xs = Krylov_Solve(Ls, bs, xs, iterParam, amgParam, amgData);
        
end

% get back the whole solution
tic
x1 = L11\(b1 - L12*xs);
toc;

% pause;
% 
% iterParam.prec_type = 'NULL';
% x = zeros(N,1); %initial guess 
% x = Krylov_Solve(L, b, x, iterParam, amgParam, amgData);
