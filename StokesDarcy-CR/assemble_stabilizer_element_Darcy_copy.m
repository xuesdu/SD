function [Sn,St] = assemble_stabilizer_element_Darcy_copy(Omega,fun,dof_trial,dof_test,P,T,neighbors,neighbors_Stokes,number_of_elements,Tb_test,Tb_trial,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    number_of_local_basis_trial,number_of_local_basis_test,basis_type_trial,basis_der_x_trial,basis_der_y_trial,index_normal_trial,...
    basis_type_test,basis_der_x_test,basis_der_y_test,index_normal_test)

Sn = sparse(dof_trial,dof_test);
St = sparse(dof_trial,dof_test);
for n = 1:number_of_elements
    if Omega == 's'
        vertices = P(:,T(:,n));
    elseif Omega == 'd'
        vertices = P(:,T(:,n+number_of_elements));
    end
    neig = neighbors_Stokes(:,n);  % 当前单元的相邻单元
    for k = 1:3
        if neig(k) ~= -1    % 这条边有相邻单元,即内部边
            if k == 1
                end_point_1 = vertices(:,2);
                end_point_2 = vertices(:,3);
            elseif k == 2
                end_point_1 = vertices(:,3);
                end_point_2 = vertices(:,1);
            elseif k == 3
                end_point_1 = vertices(:,1);
                end_point_2 = vertices(:,2);
            end
            
            neig_index = neig(k);
            if Omega == 's'
                vertices_neig = P(:,T(:,neig_index));
            elseif Omega == 'd'
                vertices_neig = P(:,T(:,neig_index+number_of_elements));
            end
            h = sqrt((end_point_1(1,1) - end_point_2(1,1))^2 + (end_point_1(2,1) - end_point_2(2,1))^2);
            mid = (end_point_1 + end_point_2)/2;
            [normal,tau] = generate_normal_vector(k,end_point_1,end_point_2);
            for alpha = 1:number_of_local_basis_trial
                for beta = 1:number_of_local_basis_test
                    % normal part
                    int_value = Gauss_quad_line_trial_test_Normal(fun,end_point_1,end_point_2,normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,vertices_neig,...
                        basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,index_normal_trial,...
                        basis_type_test,beta,basis_der_x_test,basis_der_y_test,index_normal_test);
                    i = Tb_test(beta,neig_index);
                    j = Tb_trial(alpha,n);
                    Sn(i,j) = Sn(i,j) + int_value;
                    % tangential part(midpoint formula)
                    tan_value = FE_basis_fun_local(mid(1),mid(2),vertices,basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial)*tau(index_normal_trial)*...
                        FE_basis_fun_local(mid(1),mid(2),vertices_neig,basis_type_test,beta,basis_der_x_test,basis_der_y_test)*tau(index_normal_test);
                    St(i,j) = St(i,j) + h*tan_value;
                end
            end

        else
            continue;
        end
    end
end