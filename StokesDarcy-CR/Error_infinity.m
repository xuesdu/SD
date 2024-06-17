function error = Error_infinity(solution,exact_solution,P,T,Tb_trial,Gpn,basis_type,basis_der_x,basis_der_y)

ne = size(T,2);
[Gauss_weights_ref_2D,Gauss_nodes_ref_2D] = Gauss_ref_2D(Gpn);
error = 0;
for n = 1:ne
    uh_local = solution(Tb_trial(:,n));
    vertices = P(:,T(:,n));
    [Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);

    temp = max(abs(exact_solution(Gauss_nodes_2D(:,1),Gauss_nodes_2D(:,2)) - ...
        local_FE_fun(Gauss_nodes_2D(:,1),Gauss_nodes_2D(:,2),uh_local,vertices,basis_type,basis_der_x,basis_der_y)));
    if temp > error
        error = temp;
    end
end