function [A,b] = treat_homogeneous_boundary_condition(A,b,bn_s,be_s,be_d,dof_ul)
%%%% Darcy domain: treat velocity u c\dot n = g;

global dof_Stokes

nbn = size(bn_s,2);
% P1: ul
for k = 1:nbn
    if bn_s(1,k) == -1
        i = bn_s(2,k);
        % u1
        A(i,:) = 0;
        A(:,i) = 0;
        A(i,i) = 1;
        b(i) = 0;
        % u2
        A(i+dof_ul,:) = 0;
        A(:,i+dof_ul) = 0;
        A(i+dof_ul,i+dof_ul) = 1;
        b(i+dof_ul) = 0;
    end
end

% RT0: treat uR
nbe = size(be_s,2);
for k = 1:nbe
    if be_s(6,k) == -1  % Dirichlet boundary condition
        edge = be_s(2,k);
        A(2*dof_ul+edge,:) = 0;
        A(:,2*dof_ul+edge) = 0;
        A(2*dof_ul+edge,2*dof_ul+edge) = 1;
        b(2*dof_ul+edge) = 0; 
    end
end

% RT0: treat ud
nbe = size(be_d,2);
for k = 1:nbe
    if be_d(6,k) == -1  % Dirichlet BC
        edge_local = be_d(2,k);
        if edge_local ~= 0
            A(dof_Stokes + edge_local, :) = 0;
            A(:, dof_Stokes + edge_local) = 0;
            A(dof_Stokes + edge_local, dof_Stokes +  edge_local) = 1;
            b(dof_Stokes + edge_local) =  0;
        end
    end
end




