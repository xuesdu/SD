function gb = treat_Neumann_BC_Stokes(para, P, T, be, dof_u, Tb_test, hx)


nbe = size(be,2);
e1 = sparse(dof_u,1);  e2 = sparse(dof_u,1);
for n = 1:nbe
    if be(6,n) == -2
        ele = be(3,n);
        x1 = P(1,be(4,n)); y1 = P(2,be(4,n));
        x2 = P(1,be(5,n)); y2 = P(2,be(5,n));
        x = (x1+x2)/2; y = (y1+y2)/2;
        
        if be(1,n) == 3 % top
            normal = 1;
        elseif be(1,n) == 1
            normal = -1;
        end
        
        vertices = P(:,T(:,ele));
        for alpha = 1:3  % 局部基函数的个数
            mx1 = FE_basis_fun_local(x,y,vertices,201,alpha,0,0);
            i = Tb_test(alpha,ele);
            e1(i) = e1(i) + para.gb1(x,y,para)*mx1*hx*normal;
            e2(i) = e2(i) + para.gb2(x,y,para)*mx1*hx*normal;
        end
    end
end

gb = [e1;e2];
