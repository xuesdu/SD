function out = Gauss_quad_error(uh_local,fun,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type,basis_der_x,basis_der_y)

Gpn = length(Gauss_weights_2D);
out = 0;
for n = 1:Gpn
    out = out + Gauss_weights_2D(n)*(fun(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2)) - ...
        local_FE_fun(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2),uh_local,vertices,basis_type,basis_der_x,basis_der_y))^2;
end
