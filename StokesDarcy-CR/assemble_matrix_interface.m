function A = assemble_matrix_interface(Omega,fun,dof_trial,dof_test,P,T,Inter,Tb_test,Tb_trial,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,number_of_local_basis_trial,number_of_local_basis_test,...
    basis_type_trial,basis_der_x_trial,basis_der_y_trial,...
    basis_type_test,basis_der_x_test,basis_der_y_test)

A = sparse(dof_trial,dof_test);
ele = Inter.Tf(1,:);
if Omega == 's'   % only related to the Stokes domain
    for n = 1:size(ele,2)    % loop elements
        ele_index = ele(n);  % the interface element
        vertices = P(:,T(:,ele_index));
        end_point_1 = vertices(:,3);   % two points of interface
        end_point_2 = vertices(:,1);
        if basis_type_trial == 100
            for alpha = 1:number_of_local_basis_trial
                for beta = 1:number_of_local_basis_test
                    int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,vertices,...
                        basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,...
                        basis_type_test,beta,basis_der_x_test,basis_der_y_test);
                    i = Tb_test(beta,ele_index);
                    j = Tb_trial(alpha,n);
                    A(i,j) = A(i,j) + int_value;
                end
            end
        else
            for alpha = 1:number_of_local_basis_trial
                for beta = 1:number_of_local_basis_test
                    int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,vertices,...
                        basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,...
                        basis_type_test,beta,basis_der_x_test,basis_der_y_test);
                    i = Tb_test(beta,ele_index);
                    j = Tb_trial(alpha,ele_index);
                    A(i,j) = A(i,j) + int_value;
                end
            end
        end
    end
elseif Omega == 'd' % only related to Darcy domain
    for n = 1:size(Inter.Tp,2)    
        ele_local = Inter.Tp(1,n);  % Omega_d the index starts from 1
        ele_global = Inter.Tp(2,n);  % the global index 
        vertices = P(:,T(:,ele_global));
        end_point_1 = vertices(:,3);   
        end_point_2 = vertices(:,1);
        if basis_type_trial == 100
            for alpha = 1:number_of_local_basis_trial
                for beta = 1:number_of_local_basis_test
                    int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,vertices,...
                        basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,...
                        basis_type_test,beta,basis_der_x_test,basis_der_y_test);
                    i = Tb_test(beta,ele_local);
                    j = Tb_trial(alpha,n);
                    A(i,j) = A(i,j) + int_value;
                end
            end
        else
            for alpha = 1:number_of_local_basis_trial
                for beta = 1:number_of_local_basis_test
                    int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices,vertices,...
                        basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,...
                        basis_type_test,beta,basis_der_x_test,basis_der_y_test);
                    i = Tb_test(beta,ele_local);
                    j = Tb_trial(alpha,ele_local);
                    A(i,j) = A(i,j) + int_value;
                end
            end
        end
    end
elseif Omega == 'sd'  % related to both Stokes and Darcy
    for n = 1:size(Inter.Tp,2)   
        ele_f = Inter.Tf(1,n);
        ele_p = Inter.Tp(1,n);        % Omega_p starts from 1
        ele_index_f = Inter.Tf(2,n);
        ele_index_p = Inter.Tp(2,n);  % global index
        vertices_f = P(:,T(:,ele_index_f));
        vertices_p = P(:,T(:,ele_index_p));
        end_point_1 = vertices_f(:,1);  
        end_point_2 = vertices_f(:,3);
        for alpha = 1:number_of_local_basis_trial
            for beta = 1:number_of_local_basis_test
                int_value = Gauss_quad_line_trial_test(fun,end_point_1,end_point_2,Gauss_weights_ref_1D,Gauss_nodes_ref_1D,vertices_f,vertices_p,...
                    basis_type_trial,alpha,basis_der_x_trial,basis_der_y_trial,...
                    basis_type_test,beta,basis_der_x_test,basis_der_y_test);
                i = Tb_test(beta,ele_f);
                j = Tb_trial(alpha,ele_p);
                A(i,j) = A(i,j) + int_value;
            end
        end
    end
end