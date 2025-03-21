function int_value = Gauss_quad_test_Projection(fun1,fun2,para,Gauss_weights_2D,Gauss_nodes_2D,vertices,basis_type_test,basis_index_test,basis_vector_test,basis_der_x_test,basis_der_y_test)
% 计算右端项
Gpn = length(Gauss_weights_2D);
int_value = 0;
for n = 1:Gpn
    int_value = int_value + Gauss_weights_2D(n)*dot( [fun1(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2),para),fun2(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2),para)], ...
        FE_basis_fun_local_Projection(Gauss_nodes_2D(n,1),Gauss_nodes_2D(n,2),vertices,basis_type_test,basis_index_test,basis_vector_test,basis_der_x_test,basis_der_y_test));
end