% test AMG eigensolver
% 
% @ Xiaozhe Hu, Tufts University

%---------------------
% generate Graph Laplacian
%---------------------
L = assembleGraphLaplace(100);
% A = aff;
% d = sum(A,2);
% n = size(A,1);
% D = spdiags(d,0,n,n);
% L = D - A;

%---------------------
% number of eigenvalue needed
%---------------------
number_eigen = 2;

%---------------------
% AMG parameters
%---------------------
[ amgParam ] = init_AMG_param;
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = max(10*number_eigen, 50); % size of the coarest level

amgParam.strong_connection = 0.08; 
amgParam.agg_type = 'HEC';  % 'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening
                                                    %| 'MWM': maximal weighted matching
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'
                                     
amgParam.smooth_type = 'SGS';   % 'Rd': Richardson | 'Jcb': Jacobi | 'GS': Gauss-Seidel | 'K1': Kaczmarz_AAt | 'K12': Kaczmarz_AtA
amgParam.n_presmooth = 2; % number of presmoothing      
amgParam.n_postsmooth =2; % number of postsmoothing

amgParam.tol = 1e-8;    % when AMG is used as standalone solver, tolerance for the reletive residual

% for eigenvalue problems
amgParam.number_eigen = number_eigen;

%---------------------
% AMG setup
%---------------------
amgData = AMG_Setup(L, amgParam);

% solve
tic
y = cascade_eig(amgData, 1, amgParam);
toc

%y'*L*y
eigenvalues = eigs(y'*L*y, size(y,2), 'sm');

