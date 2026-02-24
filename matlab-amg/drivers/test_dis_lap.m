% test of disaggregation
weight = 100;   % weights of the disaggregated edges
[ L_dis,  P,  Ds ] = disaggregate_laplacian( L, weight, 1 );

n = size(L,1);
n_dis = size(L_dis, 1);

Ps = Ds*P;

% from right hand side
x = rand(n,1) + (1:n)';
b = L*x;

% solve the disaggregated graph exactly
tic;
Prec_CG(L, b, zeros(n,1), @(r)disagg_exact_prec(r, L_dis, Ds, Ps), 100, 1e-6, 1);
toc;
