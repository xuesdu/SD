function A = assemble_vector_Interface_Projection(Omega,fun,dof_test,Tb_test,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_test,...
    basis_type_test,basis_vector_test,basis_der_x_test,basis_der_y_test)
global P T Inter
A = sparse(dof_test,1);
ele = Inter.Tf(1,:);
if Omega == 's'   % 只与Stokes区域有关
    for n = 1:size(ele,2)    % 遍历单元
        ele_index = ele(n);  % 界面边所在的单元
        vertices = P(:,T(:,ele_index));
        end_point_1 = vertices(:,3);   % 界面的两个端点
        end_point_2 = vertices(:,1);
        normal = [1;0];
            for beta = 1:number_of_local_basis_test
                int_value = Gauss_quad_line_trial_test_Projection(fun,end_point_1,end_point_2,normal,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,...
                    basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                i = Tb_test(beta,ele_index);
                A(i,1) = A(i,1) + int_value;
            end
    end
end