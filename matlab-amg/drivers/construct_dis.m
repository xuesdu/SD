

%n = size(A,1);
%D = spdiags(sum(A,2),0, n,n);
%L = D-A;

% get size
n = size(L,1);

% disaggregate
fprintf('----------------------------------------------------\n');
fprintf('        construct disaggregated graph \n'          );
fprintf('----------------------------------------------------\n');
tic;
[ L_dis,  P,  Ds ] = disaggregate_laplacian( L, 100, 1 );
toc;
%Ps = Ds*P;
LS_dis = Ds^(-1)*L_dis*Ds^(-1);

fprintf('----------------------------------------------------\n');
fprintf('        solve L \n'          );
fprintf('----------------------------------------------------\n');
b = ones(n,1);
b(end) = -(n-1);
testAMG(L, b);

% export
%export_MatrixMarket(L, b, './vsp_mod2_pgp2_mat.mtx',  './vsp_mod2_pgp2_rhs.mtx');

fprintf('----------------------------------------------------\n');
fprintf('        solve L_dis \n'          );
fprintf('----------------------------------------------------\n');
n_dis = size(L_dis,1);
b = ones(n_dis,1);
b(end) = -(n_dis-1);
testAMG(L_dis, b);

% export
%export_MatrixMarket(L_dis, b, './vsp_mod2_pgp2_mat_dis_w100.mtx',  './vsp_mod2_pgp2_rhs_dis_w100.mtx');

fprintf('----------------------------------------------------\n');
fprintf('        solve LS_dis \n'          );
fprintf('----------------------------------------------------\n');
b_dis = Ds^(-1)*b;
testAMG(LS_dis, b_dis);

% export
%export_MatrixMarket(LS_dis, b_dis, './vsp_mod2_pgp2_mat_scal_dis_w100.mtx',  './vsp_mod2_pgp2_rhs_scal_dis_w100.mtx');

