% Stokes-Darcy coupled problem without Lagrange multipliter
% P1+RT0: us; P0: ps; RT0: ud; P0: pd;
% 200 denotes P0 element; 201 denotes P1 element; 2011 denotes RT0 element

clear; 
global xl xr yb yt xbar
global number_of_elements number_of_edges number_of_nodes
global number_of_nodes_Stokes number_of_edges_Stokes dof_Stokes dof_Darcy
global P T E Inter alpha_T
alpha_T = 20;  % stabilizer coefficient
Gpn = 9;  

temp = 0;
for ch = 0:4   
    temp = temp+1;
    prog = select(3,ch);
    if prog.end == 1
        return
    end
    dir = prog.dir;
    para = set_parameter(dir);
    
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
    [Gauss_weights_ref_1D,Gauss_nodes_ref_1D] = Gauss_ref_1D(Gpn);
    
    [solution,iter] = assemble_matrix_vector(para,Nx,Ny,dof_ul,dof_uR,dof_ps,dof_ud,dof_pd,...
        Gauss_weights_ref_2D,Gauss_nodes_ref_2D,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    ul1 = solution(1:dof_ul);
    ul2 = solution(dof_ul+1:2*dof_ul);
    uR = solution(2*dof_ul+1:2*dof_ul+dof_uR);
    ps = solution(2*dof_ul+dof_uR+1:dof_Stokes);

    ud = solution(dof_Stokes+1:dof_Stokes+dof_ud);
    pd = solution(dof_Stokes+dof_ud+1:dof_Stokes+dof_Darcy);
    
    % L2 模
    Tb_trial_ul = T(:,1:number_of_elements/2); Eb_trial_uR = E(:,1:number_of_elements/2);
    % Stokes error
    error_L2_us1 = Error_L2_u(1,ul1,uR,para.us1,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,0,0,2011,1,0,0);
    error_L2_us2 = Error_L2_u(1,ul2,uR,para.us2,Tb_trial_ul,Eb_trial_uR,Gpn,201,0,0,0,2011,2,0,0);
    error_L2_us = sqrt(error_L2_us1^2 + error_L2_us2^2);
    error_L2_ps = Error_L2_p(1,ps,para.ps,Gpn,hx,hy);
    % Darcy error
    error_L2_ud1 = Error_L2_ud(2,para,ud,para.ud1,Eb_trial_uR,Gpn,2011,1,0,0);
    error_L2_ud2 = Error_L2_ud(2,para,ud,para.ud2,Eb_trial_uR,Gpn,2011,2,0,0);
    error_L2_ud = sqrt(error_L2_ud1^2 + error_L2_ud2^2);
    error_L2_pd = Error_L2_p(2,pd,para.pd,Gpn,hx,hy);
    
    
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