% Stokes-Darcy coupled problem without Lagrange multipliter
% P1+RT0: us; P0: ps; RT0: ud; P0: pd;
% 200 denotes P0 element; 201 denotes P1 element; 2011 denotes RT0 element

clear; 
global xl xr yb yt xbar
global number_of_elements number_of_edges number_of_nodes
global number_of_nodes_Stokes number_of_edges_Stokes dof_Stokes dof_Darcy
global P T E Inter alpha_T nu K alpha
alpha_T = 5;  % stabilizer coefficient
Gpn = 9;       % Gpn = 4 not pressure robust

temp = 0;
for ch = 1:4
    temp = temp+1;
    prog = select(1,ch);
    if prog.end == 1
        return
    end
    dir = prog.dir;
    para = set_parameter(dir);
    
    nu = para.nu(0,0); K = para.K(0,0); alpha = para.alpha(0,0);
    
    xl = para.box.left; xr = para.box.right; xbar = para.box.interface;
    yb = para.box.bottom; yt = para.box.top;
    Nx = prog.m; Ny = prog.n; hx = (xr-xl)/(2*Nx); hy = (yt-yb)/Ny;
    [P,T,E,Inter] = generate_P_T(2*Nx,Ny,201);  % 网格信息
    
    number_of_elements = size(T,2);  number_of_nodes = size(P,2); number_of_edges = (2*Nx+1)*Ny + 2*Nx*(Ny+1) + 2*Nx*Ny;
    number_of_nodes_Stokes = (Nx+1)*(Ny+1); number_of_edges_Stokes = (Nx+1)*Ny + Nx*(Ny+1) + Nx*Ny;
    
    dof_ul = number_of_nodes_Stokes; dof_uR = number_of_edges_Stokes; dof_ps = number_of_elements/2;
    dof_Stokes = dof_ul*2 + dof_uR + dof_ps;   % Stokes: total dof
    dof_ud = number_of_edges_Stokes; dof_pd = number_of_elements/2;
    dof_Darcy = dof_ud + dof_pd;   % Biot: total dof
     
    [Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
    [Gauss_weights_ref_1D,Gauss_nodes_ref_1D] = Gauss_ref_1D(4);
    
    % Assemble
    [A0, b0, A, b, Proj, solution_g] = assemble_matrix_vector_XH(para,Nx,Ny,dof_ul,dof_uR,dof_ps,dof_ud,dof_pd,...
        Gauss_weights_ref_2D,Gauss_nodes_ref_2D,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    
    % linear solve
%     iter = 0;
%     solution_0 = A\b;

    %test_stokes;
    %test_darcy;
   
    %===============================================
    % block preconditioners
    maxit = 1000;
    restart = 1000;
    tol = 1e-8;

    % reordering
    N_us = 2*dof_ul + dof_uR;
    N_ps = dof_Stokes - N_us;
    N_ud = dof_uR - Ny; 
    N_pd = dof_Darcy - Ny - N_ud; 
    Nu = N_us + N_ud;
    Np = N_ps + N_pd;

    new_dof_order = [1:N_us, dof_Stokes+1:dof_Stokes+N_ud, N_us+1:dof_Stokes, dof_Stokes+N_ud+1:dof_Stokes+N_ud+N_pd];
    [~, put_back_order] = sort(new_dof_order);
    A_reorder = A(new_dof_order, new_dof_order);
    b_reorder = b(new_dof_order);

    Auu = A_reorder(1:Nu, 1:Nu);
    Aup = A_reorder(1:Nu, Nu+1:end);
    Apu = A_reorder(Nu+1:end, 1:Nu);
    App = A_reorder(Nu+1:end, Nu+1:end);

    Mps = spdiags(0.5*hx*hy*ones(N_ps,1),0,N_ps,N_ps);
    invMps = spdiags((1/(0.5*hx*hy))*ones(N_ps,1),0,N_ps,N_ps);
    diag_invMps = (1/(0.5*hx*hy))*ones(N_ps,1);
    Mpd = spdiags(0.5*hx*hy*ones(N_pd,1),0,N_pd,N_pd);
    invMpd = spdiags((1/(0.5*hx*hy))*ones(N_pd,1),0,N_pd,N_pd);
    diag_invMpd = (1/(0.5*hx*hy))*ones(N_pd,1);

    % set parameters that can be tuned
    weight_s = 1e2; omega_s = weight_s*nu;  
    weight_d = 1e0; omega_d = weight_d*nu*K^(-1);

    Pu = Auu + Aup*[omega_s*invMps, sparse(N_ps, N_pd); sparse(N_pd, N_ps), omega_d*invMpd]*Apu;
    %Pp = [(1/omega_s)*Mps, sparse(N_ps, N_pd); sparse(N_pd, N_ps), (1/omega_d)*Mpd];
    %PA = [Pu, sparse(Nu, Np); sparse(Np, Nu), Pp]; % block diagonal
    %PA = [Pu, Aup; sparse(Np, Nu), Pp]; % block diagonal
    
    %%% ========== condition number ===================
%     A_reorder = (A_reorder+A_reorder')/2;
%     PA = (PA+PA')/2;
%     eig_max = eigs(A_reorder + (1e-12)*PA, PA, 6, 'largestabs', 'IsSymmetricDefinite', 1);
%     eig_min = eigs(A_reorder + (1e-12)*PA, PA, 6, 'smallestabs', 'IsSymmetricDefinite', 1);
%     cond_number = abs(eig_max(1))/abs(eig_min(2));

    % set AMG parameters for Pu
    amgParam = init_AMG_param;
    amgParam.print_level = 0;
    amgParam.amg_type = 'UA';
    amgParam.max_level = 2;
    amgParam.Schwarz_level = 1;

    % get blocks for Schwarz methods
    [blk_Stokes, blk_Darcy, Blk_Stokes, Blk_Darcy] = generate_block_index(Nx, Ny, dof_ul, dof_uR);
    blocks = [blk_Stokes, sparse(size(blk_Stokes,1),size(blk_Darcy,2)-(Ny+1));
              sparse(size(blk_Darcy,1),size(blk_Stokes,2)-(Ny+1)), blk_Darcy];
    %blocks = [Blk_Stokes, sparse(size(Blk_Stokes,1),size(Blk_Darcy,2)-(Ny+1));
    %         sparse(size(Blk_Darcy,1),size(Blk_Stokes,2)-(Ny+1)), Blk_Darcy];
    amgParam.Schwarz_blocks = num2cell(blocks',2);
    for i = 1:length(amgParam.Schwarz_blocks)
        [~,~,amgParam.Schwarz_blocks{i}] = find(amgParam.Schwarz_blocks{i});
    end

    % AMG setup for Pu
    amgData = AMG_Setup(Pu, amgParam);

    %[u_reorder, iter] = Prec_FGMRES(A_reorder, b_reorder, zeros(length(b),1), [], PA, maxit, restart, tol, 0);
    %[u_reorder, iter] = Prec_FGMRES(A_reorder, b_reorder, zeros(length(b),1), [], @(r)prec_diag_exact(r, Pu, diag_invMps, diag_invMpd, omega_s, omega_d), maxit, restart, tol, 0);
    [u_reorder, iter] = Prec_FGMRES(A_reorder, b_reorder, zeros(length(b),1), [], @(r)prec_diag_inexact(r, Pu, diag_invMps, diag_invMpd, omega_s, omega_d, amgParam, amgData), maxit, restart, tol, 3);
    solution_0 = u_reorder(put_back_order);
    %===============================================

    solution_0 = Proj*solution_0;
    solution = solution_0 + solution_g;

    ul1 = solution(1:dof_ul);
    ul2 = solution(dof_ul+1:2*dof_ul);
    uR = solution(2*dof_ul+1:2*dof_ul+dof_uR);
    ps = solution(2*dof_ul+dof_uR+1:dof_Stokes);

    ud = solution(dof_Stokes+1:dof_Stokes+dof_ud);
    pd = solution(dof_Stokes+dof_ud+1:dof_Stokes+dof_Darcy);
    
    fprintf("nu = %.0e, K = %0.e, cond_number = %.2f\n", nu, K, 0);
    % L2 模
    Tb_trial_ul = T(:,1:number_of_elements/2); Eb_trial_uR = E(:,1:number_of_elements/2);
    for i = 1:number_of_elements/2
        Eb_trial_p(i) = i;
    end
    % Stokes error
    error_L2_us1 = Error_L2_u(1,ul1,uR,para.us1,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,0,0,2011,1,0,0);
    error_L2_us2 = Error_L2_u(1,ul2,uR,para.us2,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,0,0,2011,2,0,0);
    error_L2_us = sqrt(error_L2_us1^2 + error_L2_us2^2);
    error_L2_ps = Error_L2_ud(1,ps,para.ps,Eb_trial_p,Gpn,200,0,0,0);
    % Darcy error
    error_L2_ud1 = Error_L2_ud(2,ud,para.ud1,Eb_trial_uR,Gpn,2011,1,0,0);
    error_L2_ud2 = Error_L2_ud(2,ud,para.ud2,Eb_trial_uR,Gpn,2011,2,0,0);
    error_L2_ud = sqrt(error_L2_ud1^2 + error_L2_ud2^2);
    error_L2_pd = K*Error_L2_ud(2,pd,para.pd,Eb_trial_p,Gpn,200,0,0,0);
    
    
    if temp == 1
        fprintf('%7.0f %s %7.2e %7.2f %s %7.2e %7.2f %s %7.2e %7.2f %s %7.2e %7.2f %s %d\n',...
            Nx,'&',error_L2_us,0,'&',error_L2_ps,0,'&',error_L2_ud,0,'&',error_L2_pd,0,'&',iter);
    else
        error_L2_us_order = log2(error_L2_us_old/error_L2_us);
        error_L2_ps_order = log2(error_L2_ps_old/error_L2_ps);
        error_L2_ud_order = log2(error_L2_ud_old/error_L2_ud);
        error_L2_pd_order = log2(error_L2_pd_old/error_L2_pd);
        fprintf('%7.0f %s %7.2e %7.2f %s %7.2e %7.2f %s %7.2e %7.2f %s %7.2e %7.2f %s %d\n',...
            Nx,'&',error_L2_us,error_L2_us_order,'&',error_L2_ps,error_L2_ps_order,...
            '&',error_L2_ud,error_L2_ud_order,'&',error_L2_pd,error_L2_pd_order,'&',iter);
    end
    error_L2_us_old = error_L2_us;
    error_L2_ps_old = error_L2_ps;
    error_L2_ud_old = error_L2_ud;
    error_L2_pd_old = error_L2_pd;

end
% draw(ul1,ul2,para,hx,hy);
% draw_p(ps,pd,para,Nx,Ny);