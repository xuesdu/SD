function ug = generate_nonhomogeneous_BC(Omega,para,bn,dof_u,dof_p,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,Nx,Ny)
global P
if Omega == 's'  % P1
    nbn_s = size(bn,2);
    ug = sparse(2*dof_u+dof_p,1);
    for k = 1:nbn_s
        if bn(2,k) ~= 0
            id = bn(2,k);
            x = P(1,id); y = P(2,id);  
            u1_exact = para.us1(x,y);
            u2_exact = para.us2(x,y);
            ug(id) = u1_exact;
            ug(id+dof_u) = u2_exact;
        end
    end  
elseif Omega == 'd'
    % Darcy boundary: ud\cdot nd = un;
    nbe_d = size(bn,2);
    ug = sparse(2*dof_u+dof_p,1);
    for k = 1:nbe_d
        if bn(1,k) ~= 0
            id = bn(2,k);
            x1 = P(1,bn(4,k)+Nx*(Ny+1)); y1 = P(2,bn(4,k)+Nx*(Ny+1));  % be(4,k) 是从头编号； +Nx*(Ny+1)为整个区域的顶点编号
            x2 = P(1,bn(5,k)+Nx*(Ny+1)); y2 = P(2,bn(5,k)+Nx*(Ny+1));
            
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
                ug(id) = exact_solution/he;
            elseif y1 == y2  % u2
                ug(id+dof_u) = exact_solution/he;
            end
        end
    end
    
end
