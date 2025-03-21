function Proj = treat_interface_condition_Strongly(para,dof_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Ny)

global Inter P T E
global dof_Stokes dof_Darcy
Proj = sparse(dof_Stokes+dof_Darcy,dof_Stokes+dof_Darcy-2*Ny);
for i = 1:dof_ul-Ny
    Proj(i,i) = 1;
end
for i = dof_ul+1:dof_Stokes
    Proj(i,i-Ny) = 1;
end

for i = dof_Stokes+Ny+1:dof_Stokes+dof_Darcy
    Proj(i,i-2*Ny) = 1;
end

for i = 1:Ny
    element_Stokes = Inter.Tf(1,i);
    element_Darcy = Inter.Tp(1,i);
    element_Darcy_global = Inter.Tp(2,i);
    
    edge_Stokes = E(2,element_Stokes);  % the number of interface deges on Stokes domain
    edge_Darcy = E(2,element_Darcy);    % the number of interface deges on Darcy domain
    
    x1 = P(1,T(1,element_Stokes));  y1 = P(2,T(1,element_Stokes));  % the endpoints' coornidates of the interface edges
    x2 = P(1,T(3,element_Stokes));  y2 = P(2,T(3,element_Stokes));
    he = abs(y1-y2);  
    
    vertices_Stokes = P(:,T(:,element_Stokes));
    vertices_Darcy = P(:,T(:,element_Darcy_global));
    
    
    % for alpha = 1:3    % P1
    %     end_point_1 = [x1;y1];   
    %     end_point_2 = [x2;y2];
    %     int_value = Gauss_quad_line_trial_test(para.one,end_point_1,end_point_2,...
    %         Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices_Stokes,...
    %         200,0,0,0,0,201,alpha,0,0,0);
    % end
    Proj(T(1,element_Stokes), dof_ul - Ny) = 1; 
    Proj(dof_Stokes+edge_Darcy,T(3,element_Stokes)) = -1;
    Proj(dof_Stokes+edge_Darcy,2*dof_ul+edge_Stokes-Ny) = 1;
end

end