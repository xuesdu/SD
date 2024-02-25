function [A0, b0, A, b, Proj, solution_g] = assemble_matrix_vector_XH(para,Nx,Ny,dof_ul,dof_uR,dof_ps,...
    dof_ud,dof_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,Gauss_weights_ref_1D,Gauss_nodes_ref_1D)

global xl xr yb yt nu K alpha
global number_of_elements 
global dof_Stokes dof_Darcy
global P T E alpha_T

hx = (xr-xl)/(2*Nx); hy = (yt-yb)/Ny;
%%%  us = ul+uR;  ud = ul; 
basis_type_trial_ul = 201;   basis_type_test_ul = 201; % P1 element
basis_type_trial_uR = 2011;  basis_type_test_uR = 2011;  % RT0 element
basis_type_trial_p = 200;    basis_type_test_p = 200; % P0 element

% ul
if basis_type_trial_ul == 201   % 有限元信息
    Pb_trial_ul = P(:,1:(Nx+1)*(Ny+1));  Tb_trial_ul = T(:,1:number_of_elements/2);
    Nx_basis_ul = Nx;   Ny_basis_ul = Ny;
    number_of_local_basis_trial_ul = 3;
end
Pb_test_ul = Pb_trial_ul;  Tb_test_ul = Tb_trial_ul; number_of_local_basis_test_ul = 3;
% uR: RT0 部分
if basis_type_trial_uR == 2011
     Eb_trial_uR = E(:,1:number_of_elements/2);
    number_of_local_basis_trial_uR = 3;
end
Eb_test_uR = Eb_trial_uR; number_of_local_basis_test_uR = 3;
% ps:trial function
Eb_trial_ps = zeros(1,number_of_elements/2);
if basis_type_trial_p == 200
    for i = 1:number_of_elements/2
        Eb_trial_ps(i) = i;  % 存储单元的信息
    end
    number_of_local_basis_trial_ps = 1;
end
Eb_test_ps = Eb_trial_ps;  number_of_local_basis_test_ps = 1;
% pd
Eb_trial_pd = Eb_trial_ps;  number_of_local_basis_trial_pd = 1; 
Eb_test_pd = Eb_trial_ps;  number_of_local_basis_test_pd = 1; 

% Omega_s
% a(ul,vl)
A1 = assemble_matrix('s',para.nu,dof_ul,dof_ul,Tb_test_ul,Tb_trial_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ul,number_of_local_basis_test_ul,basis_type_trial_ul,0,1,0,basis_type_test_ul,0,1,0);
A2 = assemble_matrix('s',para.nu,dof_ul,dof_ul,Tb_test_ul,Tb_trial_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ul,number_of_local_basis_test_ul,basis_type_trial_ul,0,0,1,basis_type_test_ul,0,0,1);
A3 = assemble_matrix('s',para.nu,dof_ul,dof_ul,Tb_test_ul,Tb_trial_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ul,number_of_local_basis_test_ul,basis_type_trial_ul,0,1,0,basis_type_test_ul,0,0,1);
% b(vl,p)
A4 = assemble_matrix('s',para.negativeone,dof_ps,dof_ul,Tb_test_ul,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ps,number_of_local_basis_test_ul,basis_type_trial_p,0,0,0,basis_type_test_ul,0,1,0);
A5 = assemble_matrix('s',para.negativeone,dof_ps,dof_ul,Tb_test_ul,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ps,number_of_local_basis_test_ul,basis_type_trial_p,0,0,0,basis_type_test_ul,0,0,1);
% a(uR,vR)
A6 = 2*assemble_matrix_Stokes_RT0('s',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,1,1,0,basis_type_test_uR,1,1,0);
A7 = 2*assemble_matrix_Stokes_RT0('s',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,2,0,1,basis_type_test_uR,2,0,1);
A6_copy = 2*assemble_matrix_Stokes_RT0('s',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,1,1,0,basis_type_test_uR,2,0,1);
A7_copy = 2*assemble_matrix_Stokes_RT0('s',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,2,0,1,basis_type_test_uR,1,1,0);
Af_RR = alpha_T*(A6+A7+A6_copy+A7_copy);
% b(vR,p)
A8 = assemble_matrix('s',para.negativeone,dof_ps,dof_uR,Eb_test_uR,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ps,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,1,1,0);
A9 = assemble_matrix('s',para.negativeone,dof_ps,dof_uR,Eb_test_uR,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ps,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,2,0,1);

% Interface terms
I1 = alpha*K^(-1/2)*assemble_matrix_interface('s',para.nu,dof_ul,dof_ul,Tb_test_ul,Tb_trial_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_ul,number_of_local_basis_test_ul,basis_type_trial_ul,1,0,0,basis_type_test_ul,1,0,0);
I2 = 2*assemble_matrix_interface('s',para.nu,dof_ul,dof_uR,Eb_test_uR,Tb_trial_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_ul,number_of_local_basis_test_uR,basis_type_trial_ul,0,1,0,basis_type_test_uR,1,0,0);

I3 = 2*assemble_matrix_interface('s',para.one,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,1,0,0,basis_type_test_uR,1,0,0);

% Omega_d 
% a(ud,vd)
Ad1 = K^(-1)*assemble_matrix('d',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,1,0,0,basis_type_test_uR,1,0,0);
Ad2 = K^(-1)*assemble_matrix('d',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,2,0,0,basis_type_test_uR,2,0,0);
% b(ud,pd)
Bd1 = assemble_matrix('d',para.negativeone,dof_pd,dof_uR,Eb_test_uR,Eb_trial_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_pd,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,1,1,0);
Bd2 = assemble_matrix('d',para.negativeone,dof_pd,dof_uR,Eb_test_uR,Eb_trial_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_pd,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,2,0,1);
% assemble total matrix
Os1 = sparse(dof_ul,dof_uR); Os2 = sparse(dof_ps,dof_ps); 
Af_LL = [2*A1+A2 A3; A3' 2*A2+A1+I1];  Af_LR = -[I2';Os1];  
% =========  gamma = -1 =================
% Aus = [Af_LL  Af_LR;
%     -Af_LR'   Af_RR];
% ========= gamma = 1; add a term gamma'(usR\cdot n, vsR\cdot n)_{Gamma}
Aus = [Af_LL  Af_LR;
      Af_LR'  Af_RR + (1e+6)*I3];
Gf_L = [A4;A5]; Gf_R = A8+A9; Bs = [Gf_L; Gf_R];
As = [Aus Bs;
     Bs' Os2];

Opd = sparse(dof_pd,dof_pd); Aud = Ad1+Ad2;  Bd = Bd1+Bd2;
Ad = [Aud Bd;
     Bd' Opd];

Asd = sparse(dof_Stokes,dof_Darcy);
A0 = [As Asd;Asd' Ad];
% assemble right hand vector
bs1 = assemble_vector_Projection('s',para.fs1,para.fs2,para,dof_ul,Tb_test_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ul,basis_type_test_ul,1,0,0);
bs2 = assemble_vector_Projection('s',para.fs1,para.fs2,para,dof_ul,Tb_test_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ul,basis_type_test_ul,2,0,0);

% bs1 = assemble_vector('s',para.fs1,para,dof_ul,Tb_test_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ul,basis_type_test_ul,0,0,0);
% bs2 = assemble_vector('s',para.fs2,para,dof_ul,Tb_test_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ul,basis_type_test_ul,0,0,0);
bs3 = assemble_vector('s',para.fs1,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,1,0,0);
bs4 = assemble_vector('s',para.fs2,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,2,0,0);
bs5 = -assemble_vector('s',para.gs,para,dof_ps,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ps,basis_type_test_p,0,0,0);

bd1 = assemble_vector('d',para.fd1,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,1,0,0);
bd2 = assemble_vector('d',para.fd2,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,2,0,0);
bd3 = -assemble_vector('d',para.gd,para,dof_pd,Eb_test_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_pd,basis_type_test_p,0,0,0);

Ib1 = assemble_vector_Interface_Projection('s',para.eta1,dof_ul,Tb_test_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_ul,basis_type_test_ul,1,0,0);
Ib2 = assemble_vector_Interface('s',para.eta1,dof_uR,Eb_test_uR,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_uR,basis_type_test_uR,1,0,0);
Ib3 = K^(-1/2)*nu*assemble_vector_Interface('s',para.eta2,dof_ul,Tb_test_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_ul,basis_type_test_ul,2,0,0);
b0 = [bs1-Ib1;bs2-Ib3;bs3+bs4-Ib2;bs5;bd1+bd2;bd3];
% treat Dirichlet boundary condition
% [bn,be] = generate_boundary_nodes_edges(2*Nx_basis_uf,Ny_basis_uf,2*Nx,Ny);
[bn_s,be_s] = generate_boundary_nodes_edges_Omegaf(Nx_basis_ul,Ny_basis_ul,Nx,Ny);
[bn_d,be_d] = generate_boundary_nodes_edges_Omegap(Nx_basis_ul,Ny_basis_ul,Nx,Ny);
%%% treat Darcy pressure p = g;
% [A,b] = treat_Dirichlet_BC(para,A,b,bn_s,be_s,be_d,Pb_trial_ul,Eb_trial_uR,dof_ul,dof_ud,Nx,Ny);
%%% treat Darcy velocity u\cdot n = g;
% solution_g = sparse(dof_Stokes+dof_Darcy,1);
% [A0,b0] = treat_boundary_condition(para,A0,b0,bn_s,be_s,be_d,Pb_trial_ul,Tb_trial_ul,dof_ul,Nx,Ny,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);

% -------------- u = u0 + ug ----------------
us1g = generate_nonhomogeneous_BC('s',para,bn_s,dof_ul,0,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny);
usRg = sparse(dof_uR,1);
udg = generate_nonhomogeneous_BC_RT0('d',para,be_d,dof_uR,0,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny);
solution_g = [us1g;usRg;sparse(dof_ps,1);udg;sparse(dof_pd,1)];
b0 = b0 - [Af_LL*us1g; Af_RR*usRg; Bs'*[us1g;usRg]; Aud*udg; Bd'*udg];
%b0 = b0 - A0*solution_g;
[A0,b0] = treat_homogeneous_boundary_condition(A0,b0,bn_s,be_s,be_d,dof_ul);

%%% treat mass conservation on the interface
Proj = treat_interface_condition(para,dof_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Ny);
A = Proj'*A0*Proj;  b = Proj'*b0;

% Additional treatment for the solution uniqueness
%A0(dof_Stokes,:) = 0;
%A0(dof_Stokes,2*dof_ul+dof_uR+1:dof_Stokes) = hx^2/2;
% A(dof_Stokes,:) = 0;
% A(dof_Stokes,2*dof_ul+dof_uR+1:dof_Stokes) = hx^2/2;
% b(dof_Stokes) = para.int_ps;

%A0(end,:) = 0;
%A0(end,dof_Stokes+dof_uR+1:end) = hx^2/2;
% A(end,:) = 0;
% A(end,dof_Stokes+(dof_uR-Ny)+1:end) = hx^2/2;
% b(end) = para.int_pd;


