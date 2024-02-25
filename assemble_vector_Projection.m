function b = assemble_vector_Projection(Omega,fun1,fun2,para,dof,Tb_test,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test,...
    basis_type_test,basis_vector_test,basis_der_x_test,basis_der_y_test)
global number_of_elements P T
b = sparse(dof,1);
if Omega == 's'
    for n = 1:number_of_elements/2 % 遍历所有单元
        vertices = P(:,T(:,n));
        [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
        for beta = 1:number_of_local_basis_test
            int_value = Gauss_quad_test_Projection(fun1,fun2,para,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
            i = Tb_test(beta,n);
            b(i) = b(i) + int_value;
        end
    end
elseif Omega == 'd'
    for n = 1:number_of_elements/2 % 遍历所有单元
        ele = n+number_of_elements/2;
        vertices = P(:,T(:,ele));
        [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
        for beta = 1:number_of_local_basis_test
            int_value = Gauss_quad_test_Projection(fun,para,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
            i = Tb_test(beta,n);
            b(i) = b(i) + int_value;
        end
    end
end