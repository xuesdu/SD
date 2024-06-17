function S = assemble_stabilizer_vector(fun,size_matrix,P,T,neighbors,number_of_elements,...
    Tb_test,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,...
    number_of_local_basis_test,basis_type_test,basis_der_x_test,basis_der_y_test)

S = sparse(size_matrix(1),1);
for n = 1:number_of_elements
    vertices = P(:,T(:,n));
    for k = 1:3  % 三角剖分（三条边）
        if neighbors(k,n) == -1
            if k == 1
                end_point_1 = vertices(:,2);
                end_point_2 = vertices(:,3);
            elseif k == 2
                end_point_1 = vertices(:,3);
                end_point_2 = vertices(:,1);
            elseif k == 3
                end_point_1 = vertices(:,1);
                end_point_2 = vertices(:,2);
            end

            for beta = 1:number_of_local_basis_test
                int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,...
                    Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,vertices,...
                    200,0,0,0,basis_type_test,beta,basis_der_x_test,basis_der_y_test);
                i = Tb_test(beta,n);
                S(i,1) = S(i,1) + int_value;
            end
        end
    end

end