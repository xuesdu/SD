function error = Error_L2(Omega,solution,fun,P,T,Tb_trial,Gpn,basis_type,basis_der_x,basis_der_y)

ne = size(T,2)/2;
[Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
error = 0;
for n = 1:ne
    uh_local = solution(Tb_trial(:,n));
    if Omega == 's'
        vertices = P(:,T(:,n));
    else
        vertices = P(:,T(:,n+ne));
    end
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);
    int_value = Gauss_quad_error(uh_local,fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type,basis_der_x,basis_der_y);
    error = error + int_value;
end

error = sqrt(error);