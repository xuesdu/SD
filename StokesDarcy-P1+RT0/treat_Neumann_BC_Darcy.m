function b = treat_Neumann_BC_Darcy(para, be, dof_u, Tb_test, Nx, Ny, hx)

global P T number_of_elements

nbe = size(be,2);

e1 = sparse(dof_u,1);  e2 = sparse(dof_u,1);
for n = 1:nbe
    if be(6,n) == -2        
        x1 = P(1,be(4,n)+Nx*(Ny+1)); y1 = P(2,be(4,n)+Nx*(Ny+1));
        x2 = P(1,be(5,n)+Nx*(Ny+1)); y2 = P(2,be(5,n)+Nx*(Ny+1));
        
        x = (x1+x2)/2; y = (y1+y2)/2;
        
        if be(1,n) == 1 % bottom boundary edge
            normal = [0;-1];
        elseif be(1,n) == 2  % right
            normal = [1;0];
        elseif be(1,n) == 3  % top
            normal = [0;1];
        else
            normal = [-1;0];
        end
        
        elements_index = be(3,n);
        vertices = P(:,T(:,elements_index+number_of_elements/2));
        for alpha = 1:3  % 局部基函数的个数
            mx1 = FE_basis_fun_local(x,y,vertices,2011,alpha,1,0,0);
            mx2 = FE_basis_fun_local(x,y,vertices,2011,alpha,2,0,0);
            i = Tb_test(alpha,elements_index);
            e1(i) = e1(i) + para.pd(x,y)*(mx1*normal(1))*hx;
            e2(i) = e2(i) + para.pd(x,y)*(mx2*normal(2))*hx;
        end
    end
end
b = e1+e2;