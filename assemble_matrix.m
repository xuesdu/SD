function A = assemble_matrix(Omega,fun,dof_trial,dof_test,Tb_test,Tb_trial,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial,number_of_local_basis_test,...
    basis_type_trial,basis_vector_trial,basis_der_x_trial,basis_der_y_trial,basis_type_test,basis_vector_test,basis_der_x_test,basis_der_y_test)

global number_of_elements P T
A = sparse(dof_test,dof_trial);
if Omega == 's'
    for n = 1:number_of_elements/2 % 遍历所有单元
        vertices = P(:,T(:,n));
        [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
        for alpha = 1:number_of_local_basis_trial
            for beta = 1:number_of_local_basis_test
                int_value = Gauss_quad_trial_test(fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,...
                    basis_type_trial,alpha,basis_vector_trial,basis_der_x_trial,basis_der_y_trial,...
                    basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                i = Tb_test(beta,n);
                j = Tb_trial(alpha,n);
                A(i,j) = A(i,j) + int_value;
            end
        end
    end
elseif Omega == 'd'
    for n = 1:number_of_elements/2 
        ele = n+number_of_elements/2;   % 单元实际的编号
        vertices = P(:,T(:,ele));
        [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
        for alpha = 1:number_of_local_basis_trial
            for beta = 1:number_of_local_basis_test
                int_value = Gauss_quad_trial_test(fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,...
                    basis_type_trial,alpha,basis_vector_trial,basis_der_x_trial,basis_der_y_trial,...
                    basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                i = Tb_test(beta,n);
                j = Tb_trial(alpha,n);
                A(i,j) = A(i,j) + int_value;
            end
        end
    end
end