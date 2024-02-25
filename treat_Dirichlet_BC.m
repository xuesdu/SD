function [A,b] = treat_Dirichlet_BC(para,A,b,bn_s,be_s,be_d,Pb_trial,Eb_trial,dof_ul,dof_ud,Nx,Ny)

global yb yt
global P T dof_Stokes
global number_of_edges_Stokes number_of_elements
hx = (yt-yb)/Nx;
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
        b(2*dof_ul+edge) = 0;
    end
end

% RT0: treat pd
e = zeros(number_of_edges_Stokes,1);
nbe = size(be_d,2);
for k = 1:nbe
    edge_local = be_d(2,k);
    if edge_local ~= 0
        x1 = P(1,be_d(4,k)+Nx*(Ny+1)); y1 = P(2,be_d(4,k)+Nx*(Ny+1));
        x2 = P(1,be_d(5,k)+Nx*(Ny+1)); y2 = P(2,be_d(5,k)+Nx*(Ny+1));
        if be_d(1,k) == 1 % bottom boundary edge
            x = (x1+x2)/2;
            y = y1;
            n1 = 0;  n2 = -1;    % normal vector n=[n1;n2];
        elseif be_d(1,k) == 2  % right
            x = x1;
            y = (y1+y2)/2;
            n1 = 1;  n2 = 0;
        elseif be_d(1,k) == 3
            x = (x1+x2)/2;
            y = y1;
            n1 = 0;  n2 = 1;
        else
            x = x1;
            y = (y1+y2)/2;
            n1 = -1;  n2 = 0;
        end
        
        
        elements_index_local = be_d(3,k);
        elements_index = elements_index_local + number_of_elements/2;
        vertices = P(:,T(:,elements_index));
        
        for alpha = 1:3
            mx1 = FE_basis_fun_local(x,y,vertices,2011,alpha,1,0,0);
            mx2 = FE_basis_fun_local(x,y,vertices,2011,alpha,2,0,0);
            i = Eb_trial(alpha,elements_index_local);
            e(i) = e(i) + (para.pd(x,y))*(mx1*n1+mx2*n2)*hx;
        end
    end
end

b(dof_Stokes + 1: dof_Stokes + dof_ud) = b(dof_Stokes + 1: dof_Stokes + dof_ud) - e;
