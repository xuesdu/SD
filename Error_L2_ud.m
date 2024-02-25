function error = Error_L2_ud(Omega,solution,fun,Eb_trial_uR,Gpn,basis_type,basis_vector,basis_der_x,basis_der_y)
global number_of_elements P T
ne = number_of_elements/2;
[Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
error = 0;
for n = 1:ne
    uh_local = solution(Eb_trial_uR(:,n));
    if Omega == 1
        vertices = P(:,T(:,n));
    elseif Omega == 2
        vertices = P(:,T(:,n+number_of_elements/2));
    end
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
    Gpn = length(Gauss_weights_2D);
    int_value = 0;
    for k = 1:Gpn
        int_value = int_value + Gauss_weights_2D(k)*(fun(Gauss_nodes_2D(k,1),Gauss_nodes_2D(k,2)) - ...
            local_FE_fun_up(Gauss_nodes_2D(k,1),Gauss_nodes_2D(k,2),uh_local,vertices,basis_type,basis_vector,basis_der_x,basis_der_y))^2;
    end
    error = error + int_value;
end

error = sqrt(error);