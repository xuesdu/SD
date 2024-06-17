function Proj = FE_basis_fun_local_Projection(x,y,vertices,basis_type,basis_index,basis_vector,basis_der_x,basis_der_y)

mid_point(:,1) = 1/2*(vertices(:,2) + vertices(:,3));
mid_point(:,2) = 1/2*(vertices(:,3) + vertices(:,1));
mid_point(:,3) = 1/2*(vertices(:,1) + vertices(:,2));
for i = 1:3
    mid_point_value(i) = FE_basis_fun_local(mid_point(1,i),mid_point(2,i),vertices,basis_type,basis_index,basis_der_x,basis_der_y);
end
normal = generate_normal_vector_Projection(vertices);
A = zeros(3);
for i = 1:3
    A(i,1) = normal(1,i);
    A(i,2) = normal(2,i);
    A(i,3) = normal(1,i)*mid_point(1,i) + normal(2,i)*mid_point(2,i);
end
basis_project_matrix(:,1) = A\[1;0;0];
basis_project_matrix(:,2) = A\[0;1;0];
basis_project_matrix(:,3) = A\[0;0;1];

q_RT1(1) = basis_project_matrix(1,1) + basis_project_matrix(3,1)*x;
q_RT1(2) = basis_project_matrix(2,1) + basis_project_matrix(3,1)*y;

q_RT2(1) = basis_project_matrix(1,2) + basis_project_matrix(3,2)*x;
q_RT2(2) = basis_project_matrix(2,2) + basis_project_matrix(3,2)*y;

q_RT3(1) = basis_project_matrix(1,3) + basis_project_matrix(3,3)*x;
q_RT3(2) = basis_project_matrix(2,3) + basis_project_matrix(3,3)*y;

p = normal(basis_vector,:).*mid_point_value;
Proj(1,1) = p(1)*q_RT1(1) + p(2)*q_RT2(1) + p(3)*q_RT3(1);
Proj(2,1) = p(1)*q_RT1(2) + p(2)*q_RT2(2) + p(3)*q_RT3(2);
end