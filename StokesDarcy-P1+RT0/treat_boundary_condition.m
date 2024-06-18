function [A,b] = treat_boundary_condition(para,A,b,bn_s,be_s,be_d,Pb_trial,Tb_trial,dof_ul,Nx,Ny,Gauss_weights_ref_1D,Gauss_nodes_ref_1D)
%%%% Darcy domain: treat velocity u c\dot n = g;

global yb yt Inter
global P T E dof_Stokes

nbn = size(bn_s,2);
% P1: ul
for k = 1:nbn
    i = bn_s(2,k);
    if i ~= 0     % boundary nodes
        % u1
        A(i,:) = 0;
        A(i,i) = 1;
        b(i) = para.us1(Pb_trial(1,i),Pb_trial(2,i));
        % u2 
        A(i+dof_ul,:) = 0;
        A(i+dof_ul,i+dof_ul) = 1;
        b(i+dof_ul) = para.us2(Pb_trial(1,i),Pb_trial(2,i));
    end
end

% RT0: treat uR
nbe = size(be_s,2);
for k = 1:nbe
    edge = be_s(2,k);
    if edge ~= 0
        A(2*dof_ul+edge,:) = 0;
        A(2*dof_ul+edge,2*dof_ul+edge) = 1;
        x1 = P(1,be_s(4,k)); y1 = P(2,be_s(4,k));
        x2 = P(1,be_s(5,k)); y2 = P(2,be_s(5,k));
        if be_s(1,k) == 2 || be_s(1,k) == 4  % vertical edge
            normal = [1;0];  he = abs(y1-y2);
            A(2*dof_ul+edge,be_s(4,k)) = 1/2;
            A(2*dof_ul+edge,be_s(5,k)) = 1/2; 
        elseif be_s(1,k) == 1 || be_s(1,k) == 3  % horizontal edge
            normal = [0;1];  he = abs(x1-x2);
            A(2*dof_ul+edge,be_s(4,k)+dof_ul) = 1/2;
            A(2*dof_ul+edge,be_s(5,k)+dof_ul) = 1/2;
        end
        exact_value = Gauss_quad_line_trial_test_normal(para.us1,para.us2,[x1;y1],[x2;y2],normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
        b(2*dof_ul+edge) = exact_value/he;
    end
end



% RT0: treat ud
nbe = size(be_d,2);
for k = 1:nbe
    edge_local = be_d(2,k);
    if edge_local ~= 0
        x1 = P(1,be_d(4,k)+Nx*(Ny+1)); y1 = P(2,be_d(4,k)+Nx*(Ny+1));
        x2 = P(1,be_d(5,k)+Nx*(Ny+1)); y2 = P(2,be_d(5,k)+Nx*(Ny+1));
        tau1 = x1 - x2;   tau2 = y1 - y2;
        if tau1 == 0  % horizontal edge
            normal = [-1;0];
            he = abs(tau2);
        elseif tau2 == 0  % vertical edge
            normal = [0;-1];
            he = abs(tau1);
        end      
        A(dof_Stokes + edge_local,:) = 0;
        A(dof_Stokes + edge_local,dof_Stokes +  edge_local) = 1;
        exact_value = Gauss_quad_line_trial_test_normal(para.ud1,para.ud2,[x1;y1],[x2;y2],normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
        b(dof_Stokes + edge_local) = exact_value/he;
    end
end




