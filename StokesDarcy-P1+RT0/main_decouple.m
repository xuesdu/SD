% Stokes-Darcy coupled problem without Lagrange multipliter
% P1+RT0: us; P0: ps; RT0: ud; P0: pd;
% 200 denotes P0 element; 201 denotes P1 element; 2011 denotes RT0 element
% generate the decoupled preconditioner 
% 3/14/2025: write by Xiaozhe

clear; 
global xl xr yb yt xbar
global number_of_elements number_of_edges number_of_nodes
global number_of_nodes_Stokes number_of_edges_Stokes dof_Stokes dof_Darcy
global P T E Inter alpha_T gamma_tilde
global nu K 
alpha_T = 20;  % stabilizer coefficient
gamma_tilde = 1e0;
Gpn = 9;  
MC = 'Weakly'; 
BC = 'I';
temp = 0;
for ch = 0:4
    temp = temp+1;
    prog = select(1,ch);
    if prog.end == 1
        return
    end
    dir = prog.dir;
    para = set_parameter(dir);
    
    xl = para.box.left; xr = para.box.right; xbar = para.box.interface;
    yb = para.box.bottom; yt = para.box.top;
    Nx = prog.m; Ny = prog.n; hx = (xr-xl)/(2*Nx); hy = (yt-yb)/Ny;
    [P,T,E,Inter] = generate_P_T(2*Nx,Ny,201);  % mesh information 

    nu = para.nu(0,0); K = para.K(0,0); alpha_BJS = para.alpha_BJS(0,0);
    if temp == 1
        fprintf("nu = %.2e, K = %.2e\n", nu, K);
    end
    
    number_of_elements = size(T,2);  number_of_nodes = size(P,2); number_of_edges = (2*Nx+1)*Ny + 2*Nx*(Ny+1) + 2*Nx*Ny;
    number_of_nodes_Stokes = (Nx+1)*(Ny+1); number_of_edges_Stokes = (Nx+1)*Ny + Nx*(Ny+1) + Nx*Ny;
    
    dof_ul = number_of_nodes_Stokes; dof_uR = number_of_edges_Stokes; dof_ps = number_of_elements/2;
    dof_Stokes = dof_ul*2 + dof_uR + dof_ps;   % Stokes: total dof
    dof_ud = number_of_edges_Stokes; dof_pd = number_of_elements/2;
    dof_Darcy = dof_ud + dof_pd;   % Biot: total dof
     
    [Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
    [Gauss_weights_ref_1D,Gauss_nodes_ref_1D] = Gauss_ref_1D(4);
    
    % Assemble
    [A0, b0, A, b, Proj, solution_g] = assemble_matrix_vector_XH(para,BC,MC,Nx,Ny,dof_ul,dof_uR,dof_ps,...
    dof_ud,dof_pd,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    %%% the nonzero elements of matrix
    % figure; spy(A); title('the non-zero elements distribution after applying the projection matrix');
    % figure; spy(A0);  title('the non-zero elements distribution');
    % figure; spy(Proj); title('Non-zero elements distribution of the preojection matrix');
    %%% linear solve
    % iter = 0;
    % solution_0 = A\b;
    % solution_0 = Proj*solution_0;
   
    maxit = 1000;
    restart = 100;
    tol = 1e-10;

    % use original ordering
    % [row_idx, col_idx, ~] = find(A(1:dof_Stokes,dof_Stokes+1:end)+A(dof_Stokes+1:end,1:dof_Stokes)');
    % u_idx = unique([row_idx;col_idx+dof_Stokes]);
    % [~, whole_idx, ~] = find(A(u_idx,:));
    % whole_idx = unique(whole_idx);
    % Att = A(whole_idx, whole_idx);
    % Nt = length(whole_idx);
    % Pt = sparse(whole_idx, (1:Nt)', ones(Nt,1), size(A,1), Nt);

    N_us = 2*dof_ul + dof_uR;
    N_ps = dof_Stokes - N_us;
    N_ud = dof_uR - Ny; 
    N_pd = dof_Darcy - Ny - N_ud; 
    Ns = N_us + N_ps;
    Nd = N_ud + N_pd;

    %%% verify the cross-term
    % C11 = A(1:N_us, Ns+1:Ns+N_ud);
    % C12 = A(1:N_us, Ns+N_ud+1:end);
    % C21 = A(N_us+1:Ns, Ns+1:Ns+N_ud);
    % C22 = A(N_us+1:Ns, Ns+N_ud+1:end);

    Mps = spdiags(0.5*hx*hy*ones(N_ps,1),0,N_ps,N_ps);
    invMps = spdiags((1/(0.5*hx*hy))*ones(N_ps,1),0,N_ps,N_ps);
    Mpd = spdiags(0.5*hx*hy*ones(N_pd,1),0,N_pd,N_pd);
    invMpd = spdiags((1/(0.5*hx*hy))*ones(N_pd,1),0,N_pd,N_pd);

    omega_s = 50*nu;
    omega_d = K^(-1);

    As_uu = A(1:N_us,1:N_us) + A(1:N_us,N_us+1:Ns)*(omega_s*invMps)*A(N_us+1:Ns,1:N_us) + A(1:N_us,Ns+N_ud+1:end)*(omega_d*invMpd)*A(Ns+N_ud+1:end,1:N_us);
    Ad_uu = A(Ns+1:Ns+N_ud,Ns+1:Ns+N_ud) + A(Ns+1:Ns+N_ud,Ns+N_ud+1:end)*(omega_d*invMpd)*A(Ns+N_ud+1:end,Ns+1:Ns+N_ud);
    Asd_u = A(1:N_us,Ns+1:Ns+N_ud) + A(1:N_us,Ns+N_ud+1:end)*(omega_d*invMpd)*A(Ns+N_ud+1:end,Ns+1:Ns+N_ud);
    Ads_u = A(Ns+1:Ns+N_ud,1:N_us) + A(Ns+1:Ns+N_ud,Ns+N_ud+1:end)*(omega_d*invMpd)*A(Ns+N_ud+1:end,1:N_us);

    %%% this is robust
    S = [As_uu,             sparse(N_us,N_ps), Asd_u,             sparse(N_us,N_pd);
         sparse(N_ps,N_us), (1/omega_s)*Mps,   sparse(N_ps,N_ud), sparse(N_ps,N_pd);
         Ads_u,             sparse(N_ud,N_ps), Ad_uu,             sparse(N_ud,N_pd);
         sparse(N_pd,N_us), sparse(N_pd,N_ps), sparse(N_pd,N_ud), (1/omega_d)*Mpd];
    %%% it is not robust
    % S = [As_uu,             sparse(N_us,N_ps), sparse(N_us,N_ud), sparse(N_us,N_pd);
    %      sparse(N_ps,N_us), (1/omega_s)*Mps,   sparse(N_ps,N_ud), sparse(N_ps,N_pd);
    %      sparse(N_ud,N_us), sparse(N_ud,N_ps), Ad_uu,             sparse(N_ud,N_pd);
    %      sparse(N_pd,N_us), sparse(N_pd,N_ps), sparse(N_pd,N_ud), (1/omega_d)*Mpd];
    
    % prec = @(r) ( Proj'*(A0\(Proj*r)) + Pt*(Att\(Pt'*r)) );
    % prec = @(r) ( Proj'*(A0\(Proj*r)) + S\r );
    prec = @(r) Proj'*(A0\(Proj*r));
    % prec = @(r) (S\r);
    %prec = @(r) asp(r, Proj, A0, S, A);
    [solution_0, iter] = Prec_FGMRES(A, b, zeros(length(b),1), [], prec, maxit, restart, tol, 0);
    solution_0 = Proj*solution_0;

    % % reordering
    % N_us = 2*dof_ul + dof_uR;
    % N_ps = dof_Stokes - N_us;
    % N_ud = dof_uR - Ny; 
    % N_pd = dof_Darcy - Ny - N_ud; 
    % Nu = N_us + N_ud;
    % Np = N_ps + N_pd;
    % 
    % new_dof_order = [1:N_us, dof_Stokes+1:dof_Stokes+N_ud, N_us+1:dof_Stokes, dof_Stokes+N_ud+1:dof_Stokes+N_ud+N_pd];
    % [~, put_back_order] = sort(new_dof_order);
    % A_reorder = A(new_dof_order, new_dof_order);
    % b_reorder = b(new_dof_order);
    % 
    % Auu = A_reorder(1:Nu, 1:Nu);
    % Aup = A_reorder(1:Nu, Nu+1:end);
    % Apu = A_reorder(Nu+1:end, 1:Nu);
    % App = A_reorder(Nu+1:end, Nu+1:end);
    % 
    % Mps = spdiags(0.5*hx*hy*ones(N_ps,1),0,N_ps,N_ps);
    % invMps = spdiags((1/(0.5*hx*hy))*ones(N_ps,1),0,N_ps,N_ps);
    % Mpd = spdiags(0.5*hx*hy*ones(N_pd,1),0,N_pd,N_pd);
    % invMpd = spdiags((1/(0.5*hx*hy))*ones(N_pd,1),0,N_pd,N_pd);
    % 
    % omega = 50*para.nu(0,0);
    % %Pu = Auu + Aup*[sparse(N_ps, N_ps), sparse(N_ps, N_pd); sparse(N_pd, N_ps), omega*invMpd]*Apu;
    % Pu = Auu + Aup*[omega*invMps, sparse(N_ps, N_pd); sparse(N_pd, N_ps), omega*invMpd]*Apu;
    % Pp = [(1/omega)*Mps, sparse(N_ps, N_pd); sparse(N_pd, N_ps), (1/omega)*Mpd];
    % PA = [Pu, sparse(Nu, Np); sparse(Np, Nu), Pp];
    % %PA = [Pu, Aup; sparse(Np, Nu), Pp];
    % %PA = [Pu, sparse(Nu, Np); Apu, Pp];
    % 
    % [u_reorder, iter] = Prec_FGMRES(A_reorder, b_reorder, zeros(length(b),1), [], PA, maxit, restart, tol, 0);
    % solution_0 = u_reorder(put_back_order);
    % solution_0 = Proj*solution_0;

    solution = solution_0 + solution_g;
    %
    ul1 = solution(1:dof_ul);
    ul2 = solution(dof_ul+1:2*dof_ul);
    uR = solution(2*dof_ul+1:2*dof_ul+dof_uR);
    ps = solution(2*dof_ul+dof_uR+1:dof_Stokes);

    ud = solution(dof_Stokes+1:dof_Stokes+dof_ud);
    pd = solution(dof_Stokes+dof_ud+1:dof_Stokes+dof_Darcy);
    
    % L2-norm
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
    error_L2_pd = Error_L2_ud(2,pd,para.pd,Eb_trial_p,Gpn,200,0,0,0);
    % Stokes H1-norm
    error_H1_us1_x = Error_L2_u(1,ul1,uR,para.us1_x,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,1,0,2011,1,1,0);
    error_H1_us1_y = Error_L2_u(1,ul1,uR,para.us1_y,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,0,1,2011,1,0,1);
    error_H1_us2_x = Error_L2_u(1,ul2,uR,para.us2_x,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,1,0,2011,2,1,0);
    error_H1_us2_y = Error_L2_u(1,ul2,uR,para.us2_y,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,0,1,2011,2,0,1);
    error_H1_us = sqrt(error_H1_us1_x^2 + error_H1_us1_y^2 + error_H1_us2_x^2 + error_H1_us2_y^2);

    if temp == 1
        fprintf('%7.0f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %d\n',...
            Nx,'&',error_H1_us,'&',0,'&',error_L2_ps,'&',0,'&',error_L2_ud,'&',0,'&',error_L2_pd,'&',0,'&',iter);        
        
    else
        error_L2_us_order = log2(error_L2_us_old/error_L2_us);
        error_H1_us_order = log2(error_H1_us_old/error_H1_us);
        error_L2_ps_order = log2(error_L2_ps_old/error_L2_ps);
        error_L2_ud_order = log2(error_L2_ud_old/error_L2_ud);
        error_L2_pd_order = log2(error_L2_pd_old/error_L2_pd);
        fprintf('%7.0f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %7.2e %s %7.2f %s %d\n',...
            Nx,'&',error_H1_us,'&',error_H1_us_order,'&',error_L2_ps,'&',error_L2_ps_order,...
            '&',error_L2_ud,'&',error_L2_ud_order,'&',error_L2_pd,'&',error_L2_pd_order,'&',iter);
    end
    error_L2_us_old = error_L2_us;
    error_H1_us_old = error_H1_us;
    error_L2_ps_old = error_L2_ps;
    error_L2_ud_old = error_L2_ud;
    error_L2_pd_old = error_L2_pd;
end

% % draw solution in Stokes domain
% P1 = P(:,1:(Nx+1)*(Ny+1));
% T1 = T(:,1:number_of_elements/2);
% draw(ul1,ul2,para,hx,hy);
% draw_p(ps,para.ps,P1,T1,number_of_elements/2,'ps');
% % draw solution in Darcy domain
% P2 = P(:,(Nx)*(Ny+1)+1:end);
% T2 = T(:,1:number_of_elements/2);
% draw_p(pd,para.pd,P2,T2,number_of_elements/2,'pd')
