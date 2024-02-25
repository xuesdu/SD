function [solution,iter] = assemble_matrix_vector(para,Nx,Ny,dof_ul,dof_uR,dof_ps,...
    dof_ud,dof_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,Gauss_weights_ref_1D,Gauss_nodes_ref_1D)

global xl xr yb yt xbar
global number_of_elements number_of_edges number_of_nodes
global number_of_nodes_Stokes number_of_edges_Stokes
global dof_Stokes dof_Darcy
global P T E Inter alpha_T

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
A6 = assemble_matrix('s',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,1,0,0,basis_type_test_uR,1,0,0);
A7 = assemble_matrix('s',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,2,0,0,basis_type_test_uR,2,0,0);
% b(vR,p)
A8 = assemble_matrix('s',para.negativeone,dof_ps,dof_uR,Eb_test_uR,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ps,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,1,1,0);
A9 = assemble_matrix('s',para.negativeone,dof_ps,dof_uR,Eb_test_uR,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial_ps,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,2,0,1);

% Interface terms
I1 = assemble_matrix_interface('s',para.one,dof_ul,dof_ul,Tb_test_ul,Tb_trial_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_ul,number_of_local_basis_test_ul,basis_type_trial_ul,1,0,0,basis_type_test_ul,1,0,0);
I2 = 2*assemble_matrix_interface('s',para.nu,dof_ul,dof_uR,Eb_test_uR,Tb_trial_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_ul,number_of_local_basis_test_uR,basis_type_trial_ul,0,1,0,basis_type_test_uR,1,0,0);

% Omega_d 
% a(ud,vd)
Ad1 = assemble_matrix('d',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,1,0,0,basis_type_test_uR,1,0,0);
Ad2 = assemble_matrix('d',para.nu,dof_uR,dof_uR,Eb_test_uR,Eb_trial_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_uR,number_of_local_basis_test_uR,basis_type_trial_uR,2,0,0,basis_type_test_uR,2,0,0);
% b(ud,pd)
Bd1 = assemble_matrix('d',para.negativeone,dof_pd,dof_uR,Eb_test_uR,Eb_trial_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_pd,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,1,1,0);
Bd2 = assemble_matrix('d',para.negativeone,dof_pd,dof_uR,Eb_test_uR,Eb_trial_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
    number_of_local_basis_trial_pd,number_of_local_basis_test_uR,basis_type_trial_p,0,0,0,basis_type_test_uR,2,0,1);
% assemble total matrix
Os1 = sparse(dof_ul,dof_uR); Os2 = sparse(dof_ps,dof_ps); 
Af_LL = [2*A1+A2 A3; A3' 2*A2+A1+I1];  Af_RR = A6+A7;  Af_LR = [I2';Os1];  
Aus = [Af_LL         Af_LR;
    -Af_LR' alpha_T*hx^(-2)*Af_RR];
Gf_L = [A4;A5]; Gf_R = A8+A9; B = [Gf_L; Gf_R];
As = [Aus B;
    -B' Os2];

Opd = sparse(dof_pd,dof_pd); Aud = Ad1+Ad2;  B = Bd1+Bd2;
Ad = [Aud B;
    -B' Opd];

Asd = sparse(dof_Stokes,dof_Darcy);
A = [As Asd;Asd' Ad];
% assemble right hand vector
bs1 = assemble_vector('s',para.fs1,para,dof_ul,Tb_test_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ul,basis_type_test_ul,0,0,0);
bs2 = assemble_vector('s',para.fs2,para,dof_ul,Tb_test_ul,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ul,basis_type_test_ul,0,0,0);
bs3 = assemble_vector('s',para.fs1,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,1,0,0);
bs4 = assemble_vector('s',para.fs2,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,2,0,0);
bs5 = assemble_vector('s',para.gs,para,dof_ps,Eb_trial_ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_ps,basis_type_test_p,0,0,0);

bd1 = assemble_vector('d',para.fd1,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,1,0,0);
bd2 = assemble_vector('d',para.fd2,para,dof_uR,Eb_test_uR,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_uR,basis_type_test_uR,2,0,0);
bd3 = assemble_vector('d',para.gd,para,dof_pd,Eb_test_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_pd,basis_type_test_p,0,0,0);

Ib1 = assemble_vector_Interface('s',para.eta1,dof_ul,Tb_test_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_ul,basis_type_test_ul,1,0,0);
Ib2 = assemble_vector_Interface('s',para.eta1,dof_uR,Eb_test_uR,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_uR,basis_type_test_uR,1,0,0);
Ib3 = assemble_vector_Interface('s',para.eta2,dof_ul,Tb_test_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_ul,basis_type_test_ul,2,0,0);
b = [bs1-Ib1;bs2-Ib3;bs3+bs4-Ib2;bs5;bd1+bd2;bd3];
% treat Dirichlet boundary condition
% [bn,be] = generate_boundary_nodes_edges(2*Nx_basis_uf,Ny_basis_uf,2*Nx,Ny);
[bn_s,be_s] = generate_boundary_nodes_edges_Omegaf(Nx_basis_ul,Ny_basis_ul,Nx,Ny);
[bn_d,be_d] = generate_boundary_nodes_edges_Omegap(Nx_basis_ul,Ny_basis_ul,Nx,Ny);
%%% treat Darcy pressure p = g;
% [A,b] = treat_Dirichlet_BC(para,A,b,bn_s,be_s,be_d,Pb_trial_ul,Eb_trial_uR,dof_ul,dof_ud,Nx,Ny);
%%% treat Darcy velocity u\cdot n = g;
[A,b] = treat_boundary_condition(para,A,b,bn_s,be_s,be_d,Pb_trial_ul,Tb_trial_ul,dof_ul,Nx,Ny,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
%%% treat mass conservation on the interface
Proj = treat_interface_condition(para,dof_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Ny);
A = Proj'*A*Proj;  b = Proj'*b;
% Additional treatment for the solution uniqueness
A(dof_Stokes,:) = 0;
A(dof_Stokes,2*dof_ul+dof_uR+1:dof_Stokes) = hx^2/2;
b(dof_Stokes) = para.int_ps;
A(end,:) = 0;
A(end,dof_Stokes+(dof_uR-Ny)+1:end) = hx^2/2;
b(end) = para.int_pd;
solution = A\b; iter = 0;
solution = Proj*solution;

