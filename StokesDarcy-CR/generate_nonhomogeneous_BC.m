function ug = generate_nonhomogeneous_BC(Omega,para,be,dof_u,dof_p,P,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny)

if Omega == 's'
    nbe_s = size(be,2);
    ug = sparse(2*dof_u+dof_p,1);
    for k = 1:nbe_s
        if be(6,k) == -1
            edge_index = be(2,k);
            x1 = P(1,be(4,k)); y1 = P(2,be(4,k));  
            x2 = P(1,be(5,k)); y2 = P(2,be(5,k));
            x_mid = (x1+x2)/2; y_mid = (y1+y2)/2;
            u1_exact = para.us1(x_mid,y_mid);
            u2_exact = para.us2(x_mid,y_mid);
            ug(edge_index) = u1_exact;
            ug(edge_index+dof_u) = u2_exact;
        end
    end  
elseif Omega == 'd'
    % Darcy boundary: ud\cdot nd = un;
    nbe_d = size(be,2);
    ug = sparse(2*dof_u+dof_p,1);
    for k = 1:nbe_d
        if be(6,k) == -1
            edge_index = be(2,k);
            x1 = P(1,be(4,k)+Nx*(Ny+1)); y1 = P(2,be(4,k)+Nx*(Ny+1));  % be(4,k) 是从头编号； +Nx*(Ny+1)为整个区域的顶点编号
            x2 = P(1,be(5,k)+Nx*(Ny+1)); y2 = P(2,be(5,k)+Nx*(Ny+1));
            
            if x1 == x2 % vertical edge
                normal = [1;0];
                he = abs(y1-y2);
            elseif y1 == y2
                normal = [0;1];
                he = abs(x1-x2);
            end
            end_point_1 = [x1;y1];  end_point_2 = [x2;y2];
            exact_solution = Gauss_quad_exact_value(para.ud1,para.ud2,end_point_1,end_point_2,normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
            if x1 == x2  % u1
                ug(edge_index) = exact_solution/he;
            elseif y1 == y2  % u2
                ug(edge_index+dof_u) = exact_solution/he;
            end
        end
    end
    
end
