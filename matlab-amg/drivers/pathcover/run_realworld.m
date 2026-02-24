% An example (real world graphs)

addpath(genpath(pwd));

% define A
%load .mat files from Stanford and UF

% get the largest connected component
G = graph(A);
S = conncomp(G);
i = mode(S);
largest_conn = find(S==i);
A = A(largest_conn,:);
A = A(:,largest_conn);
n = size(A,1);

% random initial guess for homogenous case (b=0)
rng(5)
x_init_h = rand(n,1);
% zero initial guess for nonzero rhs
x_init = zeros(n,1);

% define nonzero rhs
%low-frequency b
lfb = ones(n,1); lfb(end) = 1-n; 
%ramdom zero-sum b
rng(5)
zsb = rand(n,1);
zsb = zsb - sum(zsb)/n;

%%%%%%%%
%   Testing     %
%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard AMG test_function(A, b, tol, 'AMG', 'AMG', threshold)
%uncomment the following command to test CASE 1: b = 0
%test_function(A, zeros(n,1), x_init_h, 1e-6, 'AMG', 'AMG', 0.5)

%uncomment the following command to test CASE 2: nonzero low-frequency b
%test_function(A, lfb, x_init, 1e-6, 'AMG', 'AMG', 0.5) 

%uncomment the following command to test CASE 3: nonzero ramdom zero-sum b
%test_function(A, zsb, x_init, 1e-6, 'AMG', 'AMG', 0.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Cover Adaptive AMG test_function(A, b, tol, 'AMG_Adapt', 'NULL', threshold)
%uncomment the following command to test CASE 1: b = 0
test_function(A, zeros(n,1), x_init_h, 1e-6, 'AMG_Adapt', 'NULL',0.5)

%uncomment the following command to test CASE 2: nonzero low-frequency b
%test_function(A, lfb, x_init, 1e-6, 'AMG_Adapt', 'NULL',0.4)

%uncomment the following command to test CASE 3: nonzero ramdom zero-sum b
%test_function(A, zsb, x_init, 1e-6, 'AMG_Adapt', 'NULL',0.4)