function error = Error_L2_u(Omega,solution,uR,fun,Tb_trial,Eb_trial_uR,Gpn,basis_type,basis_vector,basis_der_x,basis_der_y,...
    basis_type_uR,basis_vector_uR,basis_der_x_uR,basis_der_y_uR)
global number_of_elements P T
ne = number_of_elements/2;
[Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
error = 0;
for n = 1:ne
    uh_local = solution(Tb_trial(:,n));
    uR_local = uR(Eb_trial_uR(:,n));
    if Omega == 1
        vertices = P(:,T(:,n));
    elseif Omega == 2
        vertices = P(:,T(:,n+number_of_elements/2));
    end
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
    int_value = Gauss_quad_error(uh_local,uR_local,fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type,basis_vector,basis_der_x,basis_der_y,...
        basis_type_uR,basis_vector_uR,basis_der_x_uR,basis_der_y_uR);
    error = error + int_value;
end

error = sqrt(error);