% nonconforming CR: Stokes-Darcy equations(CR: velocity; P0: pressure)
% 201 denotes CR element; 200 denotes P0 element
% the stabilizer only defined on the interior edges on Stokes domain
clear;
global xl xr yb yt xbar nu K alpha_BJS
Gpn = 9;    % Gauss point number
% I: no kernel; II: a kernel; III: near kernel on Darcy; IV: near kernel on Stokes;
BC = 'II';
% MC = 'Strongly';
MC = 'Weakly';  % weakly imposing the interface condition

basis_type_trial_u = 201;basis_type_test_u = 201;
basis_type_trial_p = 200;basis_type_test_p = 200;

e_us = []; e_ps = []; e_ud = []; e_pd = [];
pp = 0;
for ch = 2:6
    pp = pp+1;
    prog = select(1,ch);
    if prog.end == 1
        return
    end
    dir = prog.dir;
    para = set_parameter(dir);

    nu = para.nu(0,0);  K = para.K(0,0);  alpha_BJS = para.alpha_BJS(0,0);
    beta_t = alpha_BJS*nu^(1/2)*K^(-1/2);
    if pp == 1
        fprintf('nu = %.0e, K = %.0e \n',nu, K);
    end
    %%% ========= penalizling parameter =================
    gamma_s = nu;  gamma_d = K^(-1);  gamma_I = 100*gamma_d;

    xl = para.box.left;    xr = para.box.right;   xbar = para.box.interface;
    yb = para.box.bottom;  yt = para.box.top;
    Nx = prog.m; Ny = prog.n; hx = (xr-xl)/(2*Nx); hy = (yt-yb)/Ny; hk = sqrt(hx^2+hy^2);
    [P,T,E,neighbors,neighbors_Stokes,Inter] = generate_P_T(2*Nx,Ny);
    % single domain
    number_of_elements = size(T,2)/2;
    number_of_edges = (Nx+1)*Ny + Nx*(Ny+1) + Nx*Ny;
    % CR
    if basis_type_trial_u == 201
        Nx_basis = Nx; Ny_basis = Ny;
        Pb_trial_u = P; Tb_trial_u = T; Eb_trial_u = E;
        number_of_local_basis_trial_u = 3;
    end
    if basis_type_test_u == 201
        Nx_basis = Nx; Ny_basis = Ny;
        Pb_test_u = P; Tb_test_u = T; Eb_test_u = E;
        number_of_local_basis_test_u = 3;
    end
    Eb_trial_p = sparse(1,number_of_elements);
    % P0
    if basis_type_trial_p == 200
        for i = 1:number_of_elements
            Eb_trial_p(i) = i;
        end
        number_of_local_basis_trial_p = 1;
    end
    Eb_test_p = sparse(1,number_of_elements);
    if basis_type_test_p == 200
        for i = 1:number_of_elements
            Eb_test_p(i) = i;
        end
        number_of_local_basis_test_p = 1;
    end

    [Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
    [Gauss_weights_ref_1D,Gauss_nodes_ref_1D] = Gauss_ref_1D(4);

    %%% Dofs of unknowns
    dof_us = number_of_edges; dof_ps = number_of_elements; dof_Stokes = 2*dof_us + dof_ps;
    dof_ud = number_of_edges; dof_pd = number_of_elements; dof_Darcy = 2*dof_ud + dof_pd;

    %%% assemble matrix
    %%% Interface: <u\cdot t, v\cdot t>
    I1 = beta_t*assemble_matrix_interface('s',para.one,dof_us,dof_us,P,T,Inter,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    %%% Stokes domain
    A1 = assemble_matrix('s',para.nu,dof_us,dof_us,number_of_elements,P,T,Eb_test_u,Eb_trial_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,1,0,basis_type_test_u,1,0);
    A2 = assemble_matrix('s',para.nu,dof_us,dof_us,number_of_elements,P,T,Eb_test_u,Eb_trial_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,1,basis_type_test_u,0,1);
    A3 = assemble_matrix('s',para.nu,dof_us,dof_us,number_of_elements,P,T,Eb_test_u,Eb_trial_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,1,0,basis_type_test_u,0,1);
    % b(v,p)
    B1 = assemble_matrix('s',para.negativeone,dof_us,dof_ps,number_of_elements,P,T,Eb_test_u,Eb_trial_p,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_p,number_of_local_basis_test_u,basis_type_trial_p,0,0,basis_type_test_u,1,0);
    B2 = assemble_matrix('s',para.negativeone,dof_us,dof_ps,number_of_elements,P,T,Eb_test_u,Eb_trial_p,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_p,number_of_local_basis_test_u,basis_type_trial_p,0,0,basis_type_test_u,0,1);
    % Stabilizer: \int_e [u][v] ds; e is the interior edges
    Ss1 = assemble_stabilizer_Stokes('s',para.one,dof_us,dof_us,P,T,neighbors,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    Ss2 = assemble_stabilizer_element_Stokes('s',para.one,dof_us,dof_us,P,T,neighbors,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    Sta_s = gamma_s*[Ss1-Ss2 sparse(dof_us,dof_us);sparse(dof_us,dof_us) Ss1-Ss2]/hk;
    %%% ======= Kent: \int_e [u\cot n][v\cdot n] + [P_0(u\cdot t)][P_0(v\cdot t)] ds
    % [Sn11, St11] = assemble_stabilizer_Darcy_copy('s',para.one,dof_us,dof_us,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    %     number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,1,basis_type_test_u,0,0,1);
    % [SSn11, SSt11] = assemble_stabilizer_element_Darcy_copy('s',para.one,dof_us,dof_us,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    %     number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,1,basis_type_test_u,0,0,1);
    %
    % [Sn12,St12] = assemble_stabilizer_Darcy_copy('s',para.one,dof_us,dof_us,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    %     number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,1);
    % [SSn12,SSt12] = assemble_stabilizer_element_Darcy_copy('s',para.one,dof_us,dof_us,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    %     number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,1);
    %
    % [Sn22,St22] = assemble_stabilizer_Darcy_copy('s',para.one,dof_us,dof_us,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    %     number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,2);
    % [SSn22,SSt22] = assemble_stabilizer_element_Darcy_copy('s',para.one,dof_us,dof_us,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    %     number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,2);
    % Sn1 = gamma_s*(Sn11-SSn11)/hk; St1 = gamma_s*(St11-SSt11)/hk;
    % Sn2 = gamma_s*(Sn12-SSn12)/hk; St2 = gamma_s*(St12-SSt12)/hk;
    % Sn3 = gamma_s*(Sn22-SSn22)/hk; St3 = gamma_s*(St22-SSt22)/hk;
    % Sn = [Sn1 Sn2; Sn2' Sn3]; St = [St1 St2; St2' St3];
    % Sta_s = Sn + St;
    % =========== the jump term on the interface(on Stokes domain) ========
    Ss_I1 = assemble_matrix_interface('s',para.one,dof_us,dof_us,P,T,Inter,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    Ss_I1 = gamma_I*Ss_I1/hk;
    Aus = [2*A1+A2+Ss_I1 A3; A3' 2*A2+A1+I1];
    % Aus = [2*A1+A2    A3; A3'   2*A2+A1+I1];
    Aus = Aus + Sta_s;
    Bs = [B1;B2];  Aps = sparse(dof_ps,dof_ps);
    As = [Aus Bs; Bs' Aps];
    clear A1 A2 A3 B1 B2 Ss1 Ss2 Sta_s;

    %%% Darcy domain
    A1 = K^(-1)*assemble_matrix('d',para.one,dof_ud,dof_ud,number_of_elements,P,T,Eb_test_u,Eb_trial_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    % b(v,p)
    B1 = assemble_matrix('d',para.negativeone,dof_ud,dof_pd,number_of_elements,P,T,Eb_test_u,Eb_trial_p,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_p,number_of_local_basis_test_u,basis_type_trial_p,0,0,basis_type_test_u,1,0);
    B2 = assemble_matrix('d',para.negativeone,dof_ud,dof_pd,number_of_elements,P,T,Eb_test_u,Eb_trial_p,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,...
        number_of_local_basis_trial_p,number_of_local_basis_test_u,basis_type_trial_p,0,0,basis_type_test_u,0,1);
    % Stabilizer: sum_{E\in \Omega_d\cup \Gamma_d} int_E [u \cdot n]_E [v \cdot n]_E
    S11 = assemble_stabilizer_Darcy('d',para.one,dof_ud,dof_ud,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,1,basis_type_test_u,0,0,1);   % (u1n1)(v1n1)
    SS11 = assemble_stabilizer_element_Darcy('d',para.one,dof_ud,dof_ud,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,1,basis_type_test_u,0,0,1);
    S12 = assemble_stabilizer_Darcy('d',para.one,dof_ud,dof_ud,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,1);   % (u2n2)(v1n1)
    SS12 = assemble_stabilizer_element_Darcy('d',para.one,dof_ud,dof_ud,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,1);
    S22 = assemble_stabilizer_Darcy('d',para.one,dof_ud,dof_ud,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,2);   % (u2n2)(v2n2)
    SS22 = assemble_stabilizer_element_Darcy('d',para.one,dof_ud,dof_ud,P,T,neighbors,neighbors_Stokes,number_of_elements,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,2,basis_type_test_u,0,0,2);
    S1 = gamma_d*(S11-SS11)/hk;
    S2 = gamma_d*(S12-SS12)/hk;
    S3 = gamma_d*(S22-SS22)/hk;
    Sta = [S1 S2; S2' S3];
    % =========== the jump term on the interface(on Darcy domain) ========
    Sd_I1 = assemble_matrix_interface('d',para.one,dof_ud,dof_ud,P,T,Inter,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    Sd_I1 = gamma_I*Sd_I1/hk;
    Aud = [A1+Sd_I1 sparse(dof_ud,dof_ud);sparse(dof_ud,dof_ud) A1] + Sta;
    
    % Aud = [A1 sparse(dof_ud,dof_ud);sparse(dof_ud,dof_ud) A1] + Sta;
    Bd = [B1;B2]; Apd = sparse(dof_pd,dof_pd);
    Ad = [Aud Bd; Bd' Apd];
    clear S11 SS11 S12 SS12 S22 SS22 S1 S2 S3 A1 B1 B2;
    % =========== the jump term on the interface ========
    S_I1 = assemble_matrix_interface('sd',para.one,dof_us,dof_ud,P,T,Inter,Eb_test_u,Eb_trial_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial_u,number_of_local_basis_test_u,basis_type_trial_u,0,0,basis_type_test_u,0,0);
    S_I1 = gamma_I*S_I1/hk;
    %%% assemble right hand vector
    % Stokes domain: (f,Pi^R v)
    bs1 = assemble_vector_Projection('s',para.fs1,para.fs2,para,number_of_elements,P,T,dof_us,Eb_test_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_u,basis_type_test_u,1,0,0);
    bs2 = assemble_vector_Projection('s',para.fs1,para.fs2,para,number_of_elements,P,T,dof_us,Eb_test_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_u,basis_type_test_u,2,0,0);
    % (g,q)
    bs3 = -assemble_vector('s',para.gs,para,dof_ps,number_of_elements,P,T,Eb_test_p,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_p,basis_type_test_p,0,0);
    % interface terms
    Ib1 = assemble_vector_Interface_Projection('s',para.eta1,dof_us,P,T,Inter,Eb_test_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_u,basis_type_test_u,1,0,0);
    Ib2 = assemble_vector_Interface('s',para.eta2,dof_us,P,T,Inter,Eb_test_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
        number_of_local_basis_test_u,basis_type_test_u,0,0);
    % stabilization
    %Sb_1 = assemble_stabilizer_vector(para.us1,dof_us,P,T,neighbors,number_of_elements,Eb_test_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_test_u,basis_type_test_u,0,0);
    %Sb_2 = assemble_stabilizer_vector(para.us2,dof_us,P,T,neighbors,number_of_elements,Eb_test_u,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_test_u,basis_type_test_u,0,0);
    %bs = [bs1+gamma_s*Sb_1/hk-Ib1;bs2+gamma_s*Sb_2/hk+Ib2;bs3];
    bs = [bs1-Ib1;bs2+Ib2;bs3];
    % Darcy domain: (f,Pi^R v)
    bd1 = assemble_vector_Projection('d',para.fd1,para.fd2,para,number_of_elements,P,T,dof_ud,Eb_test_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_u,basis_type_test_u,1,0,0);
    bd2 = assemble_vector_Projection('d',para.fd1,para.fd2,para,number_of_elements,P,T,dof_ud,Eb_test_u,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_u,basis_type_test_u,2,0,0);
    bd3 = -assemble_vector('d',para.gd,para,dof_pd,number_of_elements,P,T,Eb_test_p,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test_p,basis_type_test_p,0,0);
    bd = [bd1;bd2;bd3];

    Asd = sparse(length(As),length(Ad));
    Asd(1:dof_us,1:dof_ud) = - S_I1;
    A = [ As     Asd;
        Asd'   Ad];
    %%% treat boundary condition
    [bn_s,be_s] = generate_boundary_nodes_edges_Omegas(Nx_basis,Ny_basis,Nx,Ny,BC);
    [bn_d,be_d] = generate_boundary_nodes_edges_Omegad(Nx_basis,Ny_basis,Nx,Ny,BC);

    %Stokes: Neumann
    bgs = treat_Neumann_BC_Stokes(para, P, T, be_s, dof_us, Eb_test_u, hx);
    bs = bs + [bgs; sparse(dof_ps,1)];
    %Darcy: Neumann
    bgd = treat_Neumann_BC_Darcy(para, P, T, be_d, dof_ud, Eb_test_u, number_of_elements, Nx, Ny, hx);
    bd = bd - [bgd; sparse(dof_pd,1)];
    b = [bs;bd];
    % ============== Dirichlet: us = us0 + usg; ======================
    solution_sg = generate_nonhomogeneous_BC('s',para,be_s,dof_us,dof_ps,P,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny);
    solution_dg = generate_nonhomogeneous_BC('d',para,be_d,dof_ud,dof_pd,P,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny);
    solution_g = [solution_sg;solution_dg];
    b = b - A*solution_g;
    [A,b] = treat_Stokes_Dirichlet_BC_homogeneous(A,b,be_s,P,dof_us);

    % Treat ud_0 \dot n_d = 0;
    if strcmp(BC,'II') || strcmp(BC,'III')
        Proj_S = speye(dof_Stokes);
        Proj_D = assemble_projection_matrix(neighbors, number_of_elements, be_d, Eb_trial_u, dof_ud, dof_pd, Nx, Ny);
        Proj = sparse(blkdiag(Proj_S, Proj_D));
    elseif strcmp(BC,'I') || strcmp(BC,'IV')
        Proj = speye(dof_Stokes + dof_Darcy);
    end
    A = Proj'*A*Proj; b = Proj'*b;

    % Treat mass conservation
    dof_Darcy_new = size(b,1)-dof_Stokes;  % the dofs after treating the boundary conditions
    if strcmp(MC,'Weakly')
        Proj_I = assemble_projection_matrix_interface(dof_Stokes,dof_Darcy_new,E,Inter,Ny);
    else
        Proj_I = assemble_projection_matrix_interface_Strongly(dof_Stokes,dof_Darcy_new,E,Inter,Ny);
    end
    A = Proj_I'*A*Proj_I; b = Proj_I'*b;


    %%% Treat pressure condition
    % F1: integral mean value
    A(dof_Stokes,:) = 0;  A(dof_Stokes,2*dof_us+1:dof_Stokes) = hx^2/2; b(dof_Stokes) = para.int_ps;
    A(end,:) = 0; A(end,end-dof_pd+1:end) = hx^2/2; b(end) = para.int_pd;
    % F2: fix a point value
    %exact_ps = Gauss_quad_exact_pressure(para.ps,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,[xbar xbar-hx xbar; yt yt yt-hy])/(hx^2/2);
    %exact_pd = Gauss_quad_exact_pressure(para.pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,[xr xr-hx xr; yt yt yt-hy])/(hx^2/2);
    %b_hat = [sparse(dof_Stokes-1,1);exact_ps;sparse(dof_Darcy-1,1);exact_pd];
    %b = b - A*Proj_I'*Proj'*b_hat;
    %A(dof_Stokes,:) = 0;  A(:,dof_Stokes) = 0;  A(dof_Stokes,dof_Stokes) = 1;  b(dof_Stokes) = exact_ps;
    %A(end, :) = 0;  A(:, end) = 0;  A(end, end) = 1; b(end) = exact_pd;

    % direct solve
    u0 = A\b; iter = 0; cond_number = 0;

    %     %===============================================
    %     % % Block diagnoal preconditioner
    %     maxit = 1000;  restart = 1000; tol = 1e-6;
    %     % reordering
    %     if strcmp(MC,'Weakly')
    %         N_us = 2*dof_us;
    %         N_ps = dof_Stokes - N_us;
    %         N_ud = dof_Darcy_new - dof_pd - Ny;
    %         N_pd = dof_pd;
    %         Nu = N_us + N_ud;
    %         Np = N_ps + N_pd;
    %     else
    %         N_us = 2*dof_us;
    %         N_ps = dof_Stokes - N_us;
    %         N_ud = dof_Darcy_new - dof_pd - 2*Ny;
    %         N_pd = dof_pd;
    %         Nu = N_us + N_ud;
    %         Np = N_ps + N_pd;
    %     end
    %     % new order: us ud ps pd
    %     new_dof_order = [1:N_us, dof_Stokes+1:dof_Stokes+N_ud, N_us+1:dof_Stokes, dof_Stokes+N_ud+1:dof_Stokes+N_ud+N_pd];
    %     [~,put_back_order] = sort(new_dof_order);
    %     A_reorder = A(new_dof_order,new_dof_order);
    %     b_reorder = b(new_dof_order);
    %
    %     Auu = A_reorder(1:Nu, 1:Nu);
    %     Aup = A_reorder(1:Nu, Nu+1:end);
    %     Apu = A_reorder(Nu+1:end, 1:Nu);
    %     App = A_reorder(Nu+1:end, Nu+1:end);
    %
    %     Mps = spdiags(0.5*hx*hy*ones(N_ps,1),0,N_ps,N_ps);
    %     invMps = spdiags((1/(0.5*hx*hy))*ones(N_ps,1),0,N_ps,N_ps);
    %     diag_invMps = (1/(0.5*hx*hy))*ones(N_ps,1);
    %     Mpd = spdiags(0.5*hx*hy*ones(N_pd,1),0,N_pd,N_pd);
    %     invMpd = spdiags((1/(0.5*hx*hy))*ones(N_pd,1),0,N_pd,N_pd);
    %     diag_invMpd = (1/(0.5*hx*hy))*ones(N_pd,1);
    %
    %
    %     weight_s = 100;  weight_d = 1; % a parameter that can be tuned
    %     omega_S = weight_s*nu;
    %     omega_D = weight_d*K^(-1);
    %     Pu = Auu + Aup*[omega_S*invMps, sparse(N_ps,N_pd); sparse(N_pd,N_ps), omega_D*invMpd]*Apu;
    %     Pp = [1/omega_S*Mps, sparse(N_ps,N_pd); sparse(N_pd,N_ps), 1/omega_D*Mpd];
    %     PA = [Pu, sparse(Nu,Np); sparse(Np,Nu), Pp];
    %
    %     %%% ================= Verify inf-sup condition ==================
    %     % S = Apu*(Pu)^(-1)*Aup;
    %     % I = speye(size(S));
    %     % eig_min_Schur = eigs(S + (1e-12)*I, Pp, 6, 'smallestabs');
    %     % beta_inf = sqrt(eig_min_Schur);
    %     % fprintf('beta_inf = %.2e, \n', beta_inf(1));
    %
    % %     kernel = [zeros(Nu,1); ones(N_ps,1)./2; ones(N_pd,1)];
    % %     ANS = A_reorder*kernel;
    %
    %
    %     %%% ============== condition number ===================
    %     % A_reorder = (A_reorder+A_reorder')/2;
    %     % PA = (PA+PA')/2;
    %     % eig_max = eigs(A_reorder + (1e-12)*PA, PA, 6, 'largestabs', 'IsSymmetricDefinite', 1);
    %     % eig_min = eigs(A_reorder + (1e-12)*PA, PA, 6, 'smallestabs', 'IsSymmetricDefinite', 1);
    %     % cond_number = abs(eig_max(1))/abs(eig_min(1));
    %     % cond_number_2 = abs(eig_max(1))/abs(eig_min(2));
    %     %
    %     % [V,beta_inf] = eigs(A_reorder, PA, 6, 'smallestabs', 'IsSymmetricDefinite', 1);
    %     % fprintf('cond_numer = %.2e, cond_number_eff = %.2f\n', cond_number, cond_number_2);
    %
    %     % set AMG parameters for Pu
    % %     amgParam = init_AMG_param;
    % %     amgParam.print_level = 2;
    % %     amgParam.amg_type = 'UA';
    % %     amgParam.max_level = 3;
    % %     amgParam.n_presmooth = 2; % number of presmoothing
    % %     amgParam.n_postsmooth = 2; % number of postsmoothing
    % %     amgParam.Schwarz_level = 0;
    %
    % %     % get blocks for Schwarz methods
    % %     [blk_Stokes, blk_Darcy, Blk_Stokes, Blk_Darcy] = generate_block_index(T, E, neighbors, Nx, Ny, dof_us);
    % %     blocks = [blk_Stokes, sparse(size(blk_Stokes,1),size(blk_Darcy,2)-(Ny+1));
    % %               sparse(size(blk_Darcy,1),size(blk_Stokes,2)-(Ny+1)), blk_Darcy];
    % %     %blocks = [Blk_Stokes, sparse(size(Blk_Stokes,1),size(Blk_Darcy,2)-(Ny+1));
    % %     %          sparse(size(Blk_Darcy,1),size(Blk_Stokes,2)-(Ny+1)), Blk_Darcy];
    % %     amgParam.Schwarz_blocks = num2cell(blocks',2);
    % %     for i = 1:length(amgParam.Schwarz_blocks)
    % %         [~,~,amgParam.Schwarz_blocks{i}] = find(amgParam.Schwarz_blocks{i});
    % %     end
    %
    %     % AMG setup for Pu
    %     % amgData = AMG_Setup(Pu, amgParam);
    %
    %     [u_reorder,iter,residual] = Prec_FGMRES(A_reorder, b_reorder, sparse(length(b),1), [], PA, maxit, restart, tol, 0);
    %     % [u_reorder, iter] = Prec_FGMRES(A_reorder, b_reorder, zeros(length(b),1), [], @(r)prec_diag_exact(r, Pu, diag_invMps, diag_invMpd, omega), maxit, restart, tol, 0);
    %     % [u_reorder, iter] = Prec_FGMRES(A_reorder, b_reorder, zeros(length(b),1), [], @(r)prec_diag_inexact(r, Pu, diag_invMps, diag_invMpd, omega_S, omega_D, amgParam, amgData), maxit, restart, tol, -1);
    %     u0 = u_reorder(put_back_order);
    %     % semilogy(residual,'r*--'); hold on
    %     % % ===============================================

    solution_0 = Proj*Proj_I*u0;
    solution_g = [solution_sg;solution_dg];
    solution = solution_0 + solution_g;

    us1 = solution(1:dof_us);
    us2 = solution(dof_us+1:dof_us*2);
    ps = solution(2*dof_us+1:dof_Stokes);
    ud1 = solution(dof_Stokes+1:dof_Stokes+dof_ud);
    ud2 = solution(dof_Stokes+dof_ud+1:dof_Stokes+dof_ud*2);
    pd = solution(dof_Stokes+2*dof_ud+1:dof_Stokes+dof_Darcy);

    %%% L^2-norm error
    error_L2_us1 = Error_L2('s',us1,para.us1,P,T,Eb_trial_u,Gpn,basis_type_trial_u,0,0);
    error_L2_us2 = Error_L2('s',us2,para.us2,P,T,Eb_trial_u,Gpn,basis_type_trial_u,0,0);
    error_L2_us = sqrt(error_L2_us1^2 + error_L2_us2^2);
    error_L2_ps = Error_L2('s',ps,para.ps,P,T,Eb_trial_p,Gpn,basis_type_trial_p,0,0);

    error_L2_ud1 = Error_L2('d',ud1,para.ud1,P,T,Eb_trial_u,Gpn,basis_type_trial_u,0,0);
    error_L2_ud2 = Error_L2('d',ud2,para.ud2,P,T,Eb_trial_u,Gpn,basis_type_trial_u,0,0);
    error_L2_ud = sqrt(error_L2_ud1^2 + error_L2_ud2^2);
    error_L2_pd = Error_L2('d',pd,para.pd,P,T,Eb_trial_p,Gpn,basis_type_trial_p,0,0);

    % H^1-norm
    error_H1_us1x = Error_L2('s',us1,para.us1_x,P,T,Eb_trial_u,Gpn,basis_type_trial_u,1,0);
    error_H1_us1y = Error_L2('s',us1,para.us1_y,P,T,Eb_trial_u,Gpn,basis_type_trial_u,0,1);
    error_H1_us2x = Error_L2('s',us2,para.us2_x,P,T,Eb_trial_u,Gpn,basis_type_trial_u,1,0);
    error_H1_us2y = Error_L2('s',us2,para.us2_y,P,T,Eb_trial_u,Gpn,basis_type_trial_u,0,1);
    error_H1_us = sqrt(error_H1_us1x^2 + error_H1_us1y^2 + error_H1_us2x^2 + error_H1_us2y^2);

    if pp == 1
        fprintf('%7.0f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %d\n',Nx,...
            '&',error_L2_us,'&',0,'&',error_L2_ps,'&',0,'&',error_L2_ud,'&',0,'&',error_L2_pd,'&',0,'&',iter);
    else
        error_H1_us_order = log2(error_H1_us_old/error_H1_us);
        error_L2_us_order = log2(error_L2_us_old/error_L2_us);
        error_L2_ps_order = log2(error_L2_ps_old/error_L2_ps);
        error_L2_ud_order = log2(error_L2_ud_old/error_L2_ud);
        error_L2_pd_order = log2(error_L2_pd_old/error_L2_pd);
        fprintf('%7.0f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %d\n',Nx,...
            '&',error_L2_us,'&',error_L2_us_order,'&',error_L2_ps,'&',error_L2_ps_order,...
            '&',error_L2_ud,'&',error_L2_ud_order,'&',error_L2_pd,'&',error_L2_pd_order,'&',iter);
    end

    error_L2_us_old = error_L2_us;
    error_L2_ps_old = error_L2_ps;
    error_L2_ud_old = error_L2_ud;
    error_L2_pd_old = error_L2_pd;

    error_H1_us_old = error_H1_us;

    e_us(end+1) = error_H1_us;
    e_ps(end+1) = error_L2_ps;
    e_ud(end+1) = error_L2_ud;
    e_pd(end+1) = error_L2_pd;

end
%%%

% V = Proj*Proj_I*V;
% us1 = V(1:dof_us, 1);
% us2 = V(dof_us+1:dof_us*2, 1);
% ud1 = V(dof_us*2+1:dof_us*2+dof_ud, 1);
% ud2 = V(dof_us*2+dof_ud+1:dof_us*2+dof_ud*2,1);
% ps = V(dof_us*2+dof_ud*2+1:dof_us*2+dof_ud*2+dof_ps,1);
% pd = V(dof_us*2+dof_ud*2+dof_ps+1:end,1);


% P1 = P(:,1:(Nx+1)*(Ny+1));
% draw_u(us1,para.us1,para.us1_vec,P1,T,E,number_of_elements,'us1');
% draw_u(us2,para.us2,para.us2_vec,P1,T,E,number_of_elements,'us2');
% draw_p(ps,para.ps,P1,T,number_of_elements,'ps');
% 
% P2 = P(:,(Nx)*(Ny+1)+1:end);
% draw_u(ud1,para.ud1,para.ud1_vec,P2,T,E,number_of_elements,'ud1');
% draw_u(ud2,para.ud2,para.ud1_vec,P2,T,E,number_of_elements,'ud2');
% draw_p(pd,para.pd,P2,T,number_of_elements,'pd')

fid=fopen('result_copy.txt','at');
fprintf(fid, '%.2e\t %.2e\t %.2e\t %.2e\t \n', e_us, e_ps, e_ud, e_pd);
fclose(fid);

