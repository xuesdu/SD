% test of disaggregation
weight = 10000;   % weights of the disaggregated edges
[ A_dis, P, Ds ] = disaggregate_adj( A, weight, 1 );  

n = size(A,1);
n_dis = size(A_dis, 1);

d = A*ones(n,1);
D = spdiags(d,0,n,n);
L = D - A;
d_dis = A_dis*ones(n_dis,1);
D_dis = spdiags(d_dis,0,n_dis,n_dis);
L_dis = D_dis - A_dis;

Ps = Ds*P;

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = 10; % size of the coarest level

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
amgParam.max_it = 10;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-3;    % when AMG is used as standalone solver, tolerance for the reletive residual

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother
% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'ilutp';  % nofill, crout, ilutp
amgParam.droptol = 0.01; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

%---------------------
% setup phase
%---------------------
amgData = AMG_Setup(L_dis, amgParam);

x = rand(n,1) + (1:n)';
b = L*x;

%---------------------
% solve phase
%---------------------
% solve the disaggregated graph inexactly
tic;
Prec_CG(L, b, zeros(n,1), @(r)disagg_inexact_prec(r, L_dis, Ds, Ps, amgData, amgParam), 100, 1e-6, 1);
toc;
