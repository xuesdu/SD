function b = assemble_vector(Omega,fun,para,dof,ne,P,T,Tb_test,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_test,basis_type_test,basis_der_x_test,basis_der_y_test)

b = sparse(dof,1);

for n = 1:ne % 遍历所有单元
    if Omega == 's'
        vertices = P(:,T(:,n));
    elseif Omega == 'd'
        vertices = P(:,T(:,n+ne));
    end
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
    for beta = 1:number_of_local_basis_test
        int_value = Gauss_quad_test(fun,para,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type_test,beta,basis_der_x_test,basis_der_y_test);
        i = Tb_test(beta,n);
        b(i) = b(i) + int_value;
    end
end