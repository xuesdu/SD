function ug = generate_nonhomogeneous_BC_RT0(Omega,para,be,dof_u,dof_p,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny)

global P

if Omega == 's'
    nbe_s = size(be,2);
    ug = sparse(dof_u+dof_p,1);
    for k = 1:nbe_s
        edge_index = be(2,k);
        if edge_index ~= 0
            x1 = P(1,be(4,k)); y1 = P(2,be(4,k));
            x2 = P(1,be(5,k)); y2 = P(2,be(5,k));
            if x1 == x2 % vertical edge
                normal = [1;0];
                he = abs(y1-y2);
            elseif y1 == y2
                normal = [0;1];
                he = abs(x1-x2);
            end
            end_point_1 = [x1;y1];  end_point_2 = [x2;y2];
            exact_solution = Gauss_quad_line_trial_test_normal(para.us1,para.us2,end_point_1,end_point_2,normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
            ug(edge_index) = exact_solution/he;
        end
    end  
elseif Omega == 'd'
    % Darcy boundary: ud\cdot nd = un;
    nbe_d = size(be,2);
    ug = sparse(dof_u+dof_p,1);
    for k = 1:nbe_d
        if be(1,k) ~= 0
            edge_index = be(2,k);
            x1 = P(1,be(4,k)+Nx*(Ny+1)); y1 = P(2,be(4,k)+Nx*(Ny+1));  % be(4,k) 是从头编号； +Nx*(Ny+1)为整个区域的顶点编号
            x2 = P(1,be(5,k)+Nx*(Ny+1)); y2 = P(2,be(5,k)+Nx*(Ny+1));
            
            if x1 == x2 % vertical edge
                normal = [-1;0];
                he = abs(y1-y2);
            elseif y1 == y2
                normal = [0;-1];
                he = abs(x1-x2);
            end
            end_point_1 = [x1;y1];  end_point_2 = [x2;y2];
            exact_solution = Gauss_quad_line_trial_test_normal(para.ud1,para.ud2,end_point_1,end_point_2,normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
            ug(edge_index) = exact_solution/he;
        end
    end
    
end
