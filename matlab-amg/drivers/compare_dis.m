% compare of disaggregation

%---------------------
% Generate right hand side
%---------------------
% E-R random graphs
% N = 1000;
% p = 2*log(N)/N;
% [Adj, G] = Erdos_Reyni_Random_Graph( N, p );
% L = laplacian(G); 

%A = generate_adj( '/Users/xhu03/Work/Projects/MG/Code/H-AMG/graph_examples/star-mixtures/vsp_data_and_seymourl.graph' );
%A = generate_adj( '/Users/xhu03/Work/Projects/MG/Code/H-AMG/graph_examples/star-mixtures/vsp_model1_crew1_cr42_south31.graph' );

n = size(A,1);
D = spdiags(sum(A,2),0, n,n);
L = D-A;

% get size
n = size(L,1);

% make it SPD
%[~, col] = find(L(1,:));
%Lspd = L + sparse([1;col(2);1;col(2)], [1;1;col(2);col(2)], ones(4,1), n, n);

% right hand side
%x = rand(n,1) + (1:n)';
%b = L*x;
b = ones(n,1);
b(end) = -(n-1);

% initial guess
x0 = zeros(n,1);

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
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
                                     
amgParam.smooth_type = 'GS';   % 'Rd': Richardson | 'Jb': Jacobi | 'GS': Gauss-Seidel
amgParam.n_presmooth = 1; % number of presmoothing      
amgParam.n_postsmooth =1; % number of postsmoothing
amgParam.smooth_omega = 1.2; % weight for Richardson smoother;

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 1;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-1;    % when AMG is used as standalone solver, tolerance for the reletive residual

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother
% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'ilutp';  % nofill, crout, ilutp
amgParam.droptol = 0.01; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

%---------------------
% Krylov parameters
%---------------------
max_it = 300;
tol = 1e-6; 
print_level = 1;

%---------------------
% Solve orginal problem 
%---------------------
% setup phase
amgData = AMG_Setup(L, amgParam);

% solve phase
tic
x_orig = Prec_CG(L, b, x0, @(r)AMG_prec(r, 1, amgData, amgParam), max_it, tol, print_level);
toc

%pause;

%---------------------
% Solve by disaggregted graph 
%---------------------
% disaggregate
fprintf('----------------------------------------------------\n');
fprintf('        construct disaggregated graph \n'          );
fprintf('----------------------------------------------------\n');
tic;
[ L_dis,  P,  Ds ] = disaggregate_laplacian( L, 100, 1 );
toc;
Ps = Ds*P;
%Ps = Ds*Ps;

% setup phase
%amgData_dis = AMG_Setup(L_dis, amgParam);
LS_dis = Ds^(-1)*L_dis*Ds^(-1);
amgData_dis = AMG_Setup(LS_dis, amgParam);

% solve phase
tic
x_dis = Prec_CG(L, b, x0, @(r)disagg_inexact_prec(r, L_dis, Ds, Ps, amgData_dis, amgParam), max_it, tol, print_level);
%x_dis = Prec_CG(L, b, x0, @(r)disagg_inexact_scale_prec(r, LS_dis, Ds, Ps, amgData_dis, amgParam), max_it, tol, print_level);
%x_dis = Prec_CG(L, b, x0, @(r)disagg_inexact_scale_prec(r, L_dis, Ps, amgData_dis, amgParam), max_it, tol, print_level);
toc





