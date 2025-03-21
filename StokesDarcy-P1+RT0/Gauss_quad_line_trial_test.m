function int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,...
    basis_type_trial,basis_index_trial,basis_vector_trial,basis_der_x_trial,basis_der_y_trial,...
    basis_type_test,basis_index_test,basis_vector_test,basis_der_x_test,basis_der_y_test)

Gpn = length(Gauss_weights_ref_1D);
int_value = 0;

if end_point_1(2) == end_point_2(2)  % 水平边界
    lb = min(end_point_1(1),end_point_2(1));
    ub = max(end_point_1(1),end_point_2(1));
    [Gauss_weights_1D,Gauss_nodes_1D] = Gauss_local_1D(lb,ub,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    for i = 1:Gpn
        x = Gauss_nodes_1D(i); y = end_point_1(2);
        int_value = int_value + Gauss_weights_1D(i)*fun(x,y)*...
            FE_basis_fun_local(x,y,vertices,basis_type_trial,basis_index_trial,basis_vector_trial,basis_der_x_trial,basis_der_y_trial)*...
            FE_basis_fun_local(x,y,vertices,basis_type_test,basis_index_test,basis_vector_test,basis_der_x_test,basis_der_y_test);
    end
elseif end_point_1(1) == end_point_2(1)  % 竖直边
    lb = min(end_point_1(2),end_point_2(2));
    ub = max(end_point_1(2),end_point_2(2));
    [Gauss_weights_1D,Gauss_nodes_1D] = Gauss_local_1D(lb,ub,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    for i = 1:Gpn
        x = end_point_1(1); y = Gauss_nodes_1D(i);
        int_value = int_value + Gauss_weights_1D(i)*fun(x,y)*...
            FE_basis_fun_local(x,y,vertices,basis_type_trial,basis_index_trial,basis_vector_trial,basis_der_x_trial,basis_der_y_trial)*...
            FE_basis_fun_local(x,y,vertices,basis_type_test,basis_index_test,basis_vector_test,basis_der_x_test,basis_der_y_test);
    end
else  % 斜边
    lb = min(end_point_1(1),end_point_2(1));
    ub = max(end_point_1(1),end_point_2(1));
    [Gauss_weights_1D,Gauss_nodes_1D] = Gauss_local_1D(lb,ub,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    slope = (end_point_2(2)-end_point_1(2))/(end_point_2(1)-end_point_1(1)); % 斜率
    J = sqrt(1+slope^2);
    for i = 1:Gpn
        x = Gauss_nodes_1D(i); y = slope*(x-end_point_1(1))+end_point_1(2);
        int_value = int_value + Gauss_weights_1D(i)*J*fun(x,y)*...
            FE_basis_fun_local(x,y,vertices,basis_type_trial,basis_index_trial,basis_vector_trial,basis_der_x_trial,basis_der_y_trial)*...
            FE_basis_fun_local(x,y,vertices,basis_type_test,basis_index_test,basis_vector_test,basis_der_x_test,basis_der_y_test);
    end
end


