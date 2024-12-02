function Proj = assemble_projection_matrix_interface_Strongly(dof_Stokes,dof_Darcy_new,E,Inter,Ny)

Proj = sparse(dof_Stokes+dof_Darcy_new,dof_Stokes+dof_Darcy_new-2*Ny);
for i = 1:dof_Stokes
    Proj(i,i) = 1;
end
for i = dof_Stokes+Ny+2*Ny+1+1:dof_Stokes+dof_Darcy_new
    Proj(i,i-2*Ny) = 1;
end
% special part: near the interface
for i = 1:Ny+1
    Proj(dof_Stokes+Ny+2*i-1, dof_Stokes+i) = 1;
end

% Strongly part
for i=1:Ny
    element_Stokes = Inter.Tf(1,i);
    element_Darcy = Inter.Tp(1,i);
%     element_Darcy_global = Inter.Tp(2,i);
    
    edge_Stokes_2 = E(2,element_Stokes);
    edge_Darcy_2 = E(2,element_Darcy);

    edge_Stokes_1 = E(1,element_Stokes);
    edge_Stokes_3 = E(3,element_Stokes);
    edge_Darcy_1 = E(1,element_Darcy);
    edge_Darcy_3 = E(3,element_Darcy); 

%     vertices_Stokes = P(:,T(:,element_Stokes));
%     vertices_Darcy = P(:,T(:,element_Darcy_global));
%     x1 = vertices_Stokes(1,1); y1 = vertices_Stokes(2,1);
%     x2 = vertices_Stokes(1,3); y2 = vertices_Stokes(2,3);
%     he = abs(y1-y2);
    
    Proj(dof_Stokes+edge_Darcy_2, edge_Stokes_2) = 1;   % 相当于alpha_2(只固定这一个条件即可)
    
    Proj(dof_Stokes+edge_Darcy_1, edge_Stokes_1) = -1;
    Proj(dof_Stokes+edge_Darcy_1, edge_Stokes_3) = 1;
    Proj(dof_Stokes+edge_Darcy_1, dof_Stokes+edge_Darcy_3) = 1;
    
    % double check if it is imposed strongly
%     edge_Stokes_1 = E(1,element_Stokes);
%     edge_Stokes_3 = E(3,element_Stokes);
%     edge_Darcy_1 = E(1,element_Darcy);
%     edge_Darcy_3 = E(3,element_Darcy);
%     
%     Proj(dof_Stokes + edge_Darcy_2, edge_Stokes_1) = 1; 
%     Proj(dof_Stokes + edge_Darcy_2, edge_Stokes_3) = 1;
%     Proj(dof_Stokes + edge_Darcy_2, dof_Stokes + edge_Darcy_1) = 2; 
%     Proj(dof_Stokes + edge_Darcy_2, dof_Stokes + edge_Darcy_3) = 2; 
    
end
