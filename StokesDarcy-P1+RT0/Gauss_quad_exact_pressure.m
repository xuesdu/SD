function int_value = Gauss_quad_exact_pressure(fun,Gauss_weights_ref_2D,Gauss_nodes_ref_2D,vertices)

[Gauss_weights_2D,Gauss_nodes_2D] = Gauss_local_2D(vertices,Gauss_weights_ref_2D,Gauss_nodes_ref_2D);

Gpn = length(Gauss_weights_2D);

int_value = 0;
for n = 1:Gpn
    int_value = int_value + Gauss_weights_2D(n)*fun(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2));
end