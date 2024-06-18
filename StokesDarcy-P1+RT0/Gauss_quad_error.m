function out = Gauss_quad_error(uh_local,uR_local,fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type,basis_vector,basis_der_x,basis_der_y,...
    basis_type_uR,basis_vector_uR,basis_der_x_uR,basis_der_y_uR)

Gpn = length(Gauss_weights_2D);
out = 0;
for n = 1:Gpn
    out = out + Gauss_weights_2D(n)*(fun(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2)) - ...
        local_FE_fun(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2),uh_local,uR_local,vertices,basis_type,basis_vector,basis_der_x,basis_der_y,...
        basis_type_uR,basis_vector_uR,basis_der_x_uR,basis_der_y_uR))^2;
end
