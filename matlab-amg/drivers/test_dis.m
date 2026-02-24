% test of disaggregation
weight = 100;   % weights of the disaggregated edges
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

x = rand(n,1) + (1:n)';
b = L*x;

% solve the disaggregated graph exactly
tic;
Prec_CG(L, b, zeros(n,1), @(r)disagg_exact_prec(r, L_dis, Ds, Ps), 100, 1e-6, 1);
toc;




