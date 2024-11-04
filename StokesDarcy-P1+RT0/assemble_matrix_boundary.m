function A = assemble_matrix_boundary(be,fun,dof_trial,dof_test,Tb_test,Tb_trial,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial,number_of_local_basis_test,...
    basis_type_trial,basis_vector_trial,basis_der_x_trial,basis_der_y_trial,...
    basis_type_test,basis_vector_test,basis_der_x_test,basis_der_y_test)
global P T
A = sparse(dof_test,dof_trial);

for k = 1:size(be,2)
    if be(6,k) == -2 
        if be(1,k) == 1 || be(1,k) == 2
            normal = -1;
        elseif be(1,k) == 3 || be(1,k) == 4
            normal = 1;
        end
        ele = be(3,k);
        vertices = P(:,T(:,ele));
        end_point_1 = P(:,be(4,k));
        end_point_2 = P(:,be(5,k));
        for alpha = 1:number_of_local_basis_trial
            for beta = 1:number_of_local_basis_test
                int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,...
                    basis_type_trial,alpha,basis_vector_trial,basis_der_x_trial,basis_der_y_trial,...
                    basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                i = Tb_test(beta,ele);
                j = Tb_trial(alpha,ele);
                A(i,j) = A(i,j) + int_value*normal;
            end
        end
    end
end