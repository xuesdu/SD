function Proj = treat_interface_condition(para,dof_ul,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Ny)

global Inter P T E
global dof_Stokes dof_Darcy
Proj = sparse(dof_Stokes+dof_Darcy,dof_Stokes+dof_Darcy-Ny);
for i = 1:dof_Stokes 
    Proj(i,i) = 1;
end
for i = dof_Stokes+Ny+1:dof_Stokes+dof_Darcy
    Proj(i,i-Ny) = 1;
end

for i = 1:Ny
    element_Stokes = Inter.Tf(1,i);
    element_Darcy = Inter.Tp(1,i);
    element_Darcy_global = Inter.Tp(2,i);
    
    edge_Stokes = E(2,element_Stokes);  % 界面边在Stokes区域的编号
    edge_Darcy = E(2,element_Darcy);    % 界面边在Darcy区域的编号
    
    x1 = P(1,T(1,element_Stokes));  y1 = P(2,T(1,element_Stokes));  % 界面边的两个端点坐标
    x2 = P(1,T(3,element_Stokes));  y2 = P(2,T(3,element_Stokes));
    he = abs(y1-y2);  
    
    vertices_Stokes = P(:,T(:,element_Stokes));
    vertices_Darcy = P(:,T(:,element_Darcy_global));
    
    
    for alpha = 1:3    % P1
        end_point_1 = [x1;y1];   % 界面的两个端点
        end_point_2 = [x2;y2];
        int_value = Gauss_quad_line_trial_test(para.one,end_point_1,end_point_2,...
            Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices_Stokes,...
            200,0,0,0,0,201,alpha,0,0,0);
        Proj(dof_Stokes+edge_Darcy,T(alpha,element_Stokes)) = -int_value/he;
        %Proj(dof_Stokes+edge_Darcy,T(alpha,element_Stokes)) = 0;  % P1部分系数取0(P1 不起作用)
    end
    Proj(dof_Stokes+edge_Darcy,2*dof_ul+edge_Stokes) = 1;
    %Proj(dof_Stokes+edge_Darcy,2*dof_ul+edge_Stokes) = -1;  % RT0 部分正好抵消
end

end