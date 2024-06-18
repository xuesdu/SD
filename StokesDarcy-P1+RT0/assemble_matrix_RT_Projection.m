function Proj = assemble_matrix_RT_Projection(para,dof_ul,dof_uR,...
    number_of_local_basis_test,basis_type_test,basis_vector_test,basis_der_x_test,basis_der_y_test)
global P T E number_of_elements

Proj = sparse(2*dof_ul,dof_uR);
for n = 1:number_of_elements/2
    vertices = P(:,T(:,n));    
    for beta = 1:number_of_local_basis_test
        id_ul = T(beta,n);
        for k = 1:3        % 遍历三角形的边
            id_uR = E(k,n);
            switch k
                case 1
                    mid_1 = (vertices(:,2) + vertices(:,3))/2;
                    tau = [vertices(1,2) - vertices(1,3); vertices(2,2) - vertices(2,3)];
                    normal = [tau(2);-tau(1)];
                    if dot(normal,[1;1]) > 0
                        normal = normal;
                    else
                        normal = -normal;
                    end
                    basis_value = FE_basis_fun_local(mid_1(1),mid_1(2),vertices,basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                    value1 = basis_value*normal(1);
                    value2 = basis_value*normal(2);
                    Proj(id_ul,id_uR) = Proj(id_ul,id_uR) +  value1;
                    Proj(id_ul+dof_ul,id_uR) = Proj(id_ul+dof_ul,id_uR) +  value2;
                case 2
                    mid_2 = (vertices(:,3) + vertices(:,1))/2;
                    tau = [vertices(1,3) - vertices(1,1); vertices(2,3) - vertices(2,1)];
                    normal = [tau(2);-tau(1)];
                    if dot(normal,[-1;0]) > 0
                        normal = normal;
                    else
                        normal = -normal;
                    end
                    basis_value = FE_basis_fun_local(mid_2(1),mid_2(2),vertices,basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                    value1 = basis_value*normal(1);
                    value2 = basis_value*normal(2);
                    Proj(id_ul,id_uR) = Proj(id_ul,id_uR) + value1;
                    Proj(id_ul+dof_ul,id_uR) = Proj(id_ul+dof_ul,id_uR) + value2;
                case 3
                    mid_3 = (vertices(:,1) + vertices(:,2))/2;
                    tau = [vertices(1,1) - vertices(1,2); vertices(2,1) - vertices(2,2)];
                    normal = [tau(2);-tau(1)];
                    if dot(normal,[0;-1]) > 0
                        normal = normal;
                    else
                        normal = -normal;
                    end
                    basis_value = FE_basis_fun_local(mid_3(1),mid_3(2),vertices,basis_type_test,beta,basis_vector_test,basis_der_x_test,basis_der_y_test);
                    value1 = basis_value*normal(1);
                    value2 = basis_value*normal(2);
                    Proj(id_ul,id_uR) = Proj(id_ul,id_uR) + value1;
                    Proj(id_ul+dof_ul,id_uR) = Proj(id_ul+dof_ul,id_uR) + value2;
            end
        end
    end

end