function error = Error_infinity_u(solution,uR,fun,P,T,Tb_trial,Eb_trial,Gpn,basis_type,basis_vector,basis_der_x,basis_der_y,...
    basis_type_uR,basis_vector_uR,basis_der_x_uR,basis_der_y_uR)

ne = size(T,2);
[Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
error = 0;
for n = 1:ne
    uh_local = solution(Tb_trial(:,n));
    uR_local = uR(Eb_trial(:,n));
    vertices = P(:,T(:,n));
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
    temp = max(abs(fun(Gauss_nodes_2D(:,1),Gauss_nodes_2D(:,2)) - ...
        (local_FE_fun(Gauss_nodes_2D(:,1),Gauss_nodes_2D(:,2),uh_local,vertices,basis_type,basis_vector,basis_der_x,basis_der_y)+...
        local_FE_fun(Gauss_nodes_2D(:,1),Gauss_nodes_2D(:,2),uR_local,vertices,basis_type_uR,basis_vector_uR,basis_der_x_uR,basis_der_y_uR))));
    if temp > error
        error = temp;
    end
end