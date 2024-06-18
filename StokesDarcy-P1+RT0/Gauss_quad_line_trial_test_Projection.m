function int_value = Gauss_quad_line_trial_test_Projection(fun,end_point_1,end_point_2,normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,...
    basis_type,basis_index,basis_vector,basis_der_x,basis_der_y)

Gpn = length(Gauss_weights_ref_1D);
int_value = 0;

if end_point_1(2) == end_point_2(2)  % 水平边界
    lb = min(end_point_1(1),end_point_2(1));
    ub = max(end_point_1(1),end_point_2(1));
    [Gauss_weights_1D,Gauss_nodes_1D] = Gauss_local_1D(lb,ub,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    for i = 1:Gpn
        x = Gauss_nodes_1D(i); y = end_point_1(2);
        int_value = int_value + Gauss_weights_1D(i)*fun(x,y)*...
            dot(FE_basis_fun_local_Projection(x,y,vertices,edge_index,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y),normal);
    end
elseif end_point_1(1) == end_point_2(1)  % 竖直边
    lb = min(end_point_1(2),end_point_2(2));
    ub = max(end_point_1(2),end_point_2(2));
    [Gauss_weights_1D,Gauss_nodes_1D] = Gauss_local_1D(lb,ub,Gauss_weights_ref_1D,Gauss_nodes_ref_1D);
    for i = 1:Gpn
        x = end_point_1(1); y = Gauss_nodes_1D(i);
        int_value = int_value + Gauss_weights_1D(i)*fun(x,y)*...
            dot(FE_basis_fun_local_Projection(x,y,vertices,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y),normal);
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
            dot(FE_basis_fun_local(x,y,vertices,edge_index,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y),normal);
    end
end


