% An example (regular grids & ring graphs & ER random graphs)

addpath(genpath(pwd));

% define A
%load('CA-GrQc.mat')
%L_1 = L;
%n_b = length(L); %multiply P_mat by b to get the extended rhs


n_l = 2^6;
%level = 2;
level = log2(n_l^2)
% regular grid
L_1 = assembleGraphLaplace(n_l); %create a regular grid
% ring graph
%[~,L_1] = WattsStrogatz(n_l^2,2,0);
%[~,L_1] = WattsStrogatz(n_l^2,2,1/n_l); %create a ring graph with
% ER random graph
%[ ~, G ] = Erdos_Reyni_Random_Graph( n_l^2, 5*log(n_l^2)/n_l^2 );
%L_1 = laplacian(G);

n_b = n_l^2; % size of L_1 - original graph



%example_SPT;
[L_T, L_M_rel, L_tilde, ind_reorder,P_mat, layer_ind] = build_MSP(L_1, 1/sqrt(n_b), level,'MWM');
H=L_T(ind_reorder,ind_reorder); % reorder index by aggregations to avoid adding fill-ins
L = L_tilde(ind_reorder,ind_reorder);
M = L_M_rel(ind_reorder,ind_reorder);
n = size(H,1);

% % random initial guess for homogenous case (b=0)
%rng(5)
%x_init_h = rand(n,1);
% b = L*x_init_h;
% zero initial guess for nonzero rhs
%x_init = zeros(n,1);
x_init = zeros(n_b,1);

% define nonzero rhs
%low-frequency b
lfb = ones(n_b,1); lfb(end) = 1-n_b; 
%lfb = P_mat*lfb;
%ramdom zero-sum b
rng(5)
zsb = rand(n_b,1);
zsb = zsb - sum(zsb)/n_b;
%zsb = P_mat*zsb;

%%%%%%%%
%   Testing     %
%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard AMG test_function(A, b, x_init, tol, solver_type, prec_type, threshold)
%uncomment the following command to test CASE 1: b = 0
%test_function(L, A, x_init_h, x_init, 1e-8, 'AMG', 'AMG', 0.5, layer_ind, level) 
%test_function(L, A, x_init_h, x_init, 1e-8, 'CG', 'MSP', 0.5) 

%uncomment the following command to test CASE 2: nonzero low-frequency b
%test_function(L, M, M, lfb, x_init, 1e-8, 'CG', 'AMG', 0.5, layer_ind, level) 

% using positive subgraph to solve original
test_function(L_1, [], M, lfb, x_init, 1e-8, 'CG', '', 0.5,layer_ind,level, P_mat) 
% using MSP on positive subgraph to solve original
test_function(L_1, H, M, lfb, x_init, 1e-8, 'FGMRES', 'MSP', 0.5, layer_ind, level, P_mat) 
% using MSP directly to solve original
test_function(L_1, H, [], lfb, x_init, 1e-8, 'CG', '', 0.5, layer_ind, level, P_mat) 

%uncomment the following command to test CASE 3: nonzero ramdom zero-sum b
%test_function(A, zsb, x_init, 1e-8, 'AMG', 'AMG', 0.5)
%test_function(A, A, zsb, x_init, 1e-8, 'CG', 'AMG', 0.5, layer_ind, level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Cover Adaptive AMG test_function(A, b, x_init, tol, solver_type, prec_type, threshold)
%uncomment the following command to test CASE 1: b = 0
%test_function(A, zeros(n,1), x_init_h, 1e-8, 'AMG_Adapt', 'NULL',0.5)

%uncomment the following command to test CASE 2: nonzero low-frequency b
%test_function(A, lfb, x_init, 1e-8, 'AMG_Adapt', 'NULL',0.4)

%uncomment the following command to test CASE 3: nonzero ramdom zero-sum b
%test_function(A, zsb, x_init, 1e-8, 'AMG_Adapt', 'NULL',0.4)