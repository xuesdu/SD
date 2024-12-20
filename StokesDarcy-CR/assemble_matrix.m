function A = assemble_matrix(Omega,fun,dof_trial,dof_test,ne,P,T,Eb_test,Eb_trial,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,number_of_local_basis_trial,number_of_local_basis_test,basis_type_trial,basis_der_x_trial,basis_der_y_trial,basis_type_test,basis_der_x_test,basis_der_y_test)

A = sparse(dof_trial,dof_test);

for n = 1:ne 
    if Omega == 's'
        vertices = P(:,T(:,n));
    elseif Omega == 'd'
        vertices = P(:,T(:,n+ne));
    end
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
    for alpha = 1:number_of_local_basis_trial
        for beta = 1:number_of_local_basis_test
            int_value = Gauss_quad_trial_test(fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,basis_type_test,beta,basis_der_x_test,basis_der_y_test);
            i = Eb_test(beta,n);
            j = Eb_trial(alpha,n);
            A(i,j) = A(i,j) + int_value;
        end
    end
end