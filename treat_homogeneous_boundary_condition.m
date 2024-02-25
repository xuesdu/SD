function [A,b] = treat_homogeneous_boundary_condition(A,b,bn_s,be_s,be_d,dof_ul)
%%%% Darcy domain: treat velocity u c\dot n = g;

global dof_Stokes

nbn = size(bn_s,2);
% P1: ul
for k = 1:nbn
    i = bn_s(2,k);
    if i ~= 0     % boundary nodes
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
    edge = be_s(2,k);
    if edge ~= 0
        A(2*dof_ul+edge, :) = 0;
        A(:, 2*dof_ul+edge) = 0;
        A(2*dof_ul+edge, 2*dof_ul+edge) = 1;
%         x1 = P(1,be_s(4,k)); y1 = P(2,be_s(4,k));
%         x2 = P(1,be_s(5,k)); y2 = P(2,be_s(5,k));
%         if be_s(1,k) == 2 || be_s(1,k) == 4  % vertical edge
%             normal = [1;0];  he = abs(y1-y2);
%             A(2*dof_ul+edge,be_s(4,k)) = 1/2;
%             A(2*dof_ul+edge,be_s(5,k)) = 1/2; 
%         elseif be_s(1,k) == 1 || be_s(1,k) == 3  % horizontal edge
%             normal = [0;1];  he = abs(x1-x2);
%             A(2*dof_ul+edge,be_s(4,k)+dof_ul) = 1/2;
%             A(2*dof_ul+edge,be_s(5,k)+dof_ul) = 1/2;
%         end
%         exact = Gauss_quad_line_trial_test_normal(para.us1,para.us2,[x1;y1],[x2;y2],normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
        b(2*dof_ul+edge) = 0;    % exact/he;
    end
end



% RT0: treat ud
nbe = size(be_d,2);
for k = 1:nbe
    edge_local = be_d(2,k);
    if edge_local ~= 0
        A(dof_Stokes + edge_local, :) = 0;
        A(:, dof_Stokes + edge_local) = 0;
        A(dof_Stokes + edge_local,dof_Stokes +  edge_local) = 1;
        b(dof_Stokes + edge_local) =  0;
    end
end




