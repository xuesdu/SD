function [e1, e2] = treat_Neumann_BC_Stokes(para, be, dof_u, Tb_test, hx, basis_type)

global P T 
nbe = size(be,2);
e1 = sparse(dof_u,1);  e2 = sparse(dof_u,1);
for n = 1:nbe
    if be(6,n) == -2
        ele = be(3,n);
        x1 = P(1,be(4,n)); y1 = P(2,be(4,n));
        x2 = P(1,be(5,n)); y2 = P(2,be(5,n));
        x = (x1+x2)/2; y = (y1+y2)/2;
        w = hx/6;
        
        if be(1,n) == 3 % top
            normal = 1;
        elseif be(1,n) == 1
            normal = -1;
        end
        
        vertices = P(:,T(:,ele));
        for alpha = 1:3  % 局部基函数的个数
            mx1 = FE_basis_fun_local(x,y,vertices,basis_type,alpha,1,0,0);
            mx2 = FE_basis_fun_local(x,y,vertices,basis_type,alpha,2,0,0);
            i = Tb_test(alpha,ele);
            e1(i) = e1(i) + para.gb1(x,y,para)*mx1*hx*normal;
            e2(i) = e2(i) + para.gb2(x,y,para)*mx2*hx*normal;
        end
        

%         for alpha = 1:3  % 局部基函数的个数
%             mx11 = FE_basis_fun_local(x1,y1,vertices,basis_type,alpha,1,0,0);
%             mx12 = FE_basis_fun_local(x,y,vertices,basis_type,alpha,1,0,0);
%             mx13 = FE_basis_fun_local(x2,y2,vertices,basis_type,alpha,1,0,0);
%             
%             mx21 = FE_basis_fun_local(x1,y1,vertices,basis_type,alpha,2,0,0);
%             mx22 = FE_basis_fun_local(x,y,vertices,basis_type,alpha,2,0,0);
%             mx23 = FE_basis_fun_local(x2,y2,vertices,basis_type,alpha,2,0,0);
%             
%             i = Tb_test(alpha,ele);
%             e1(i) = e1(i) + w*(para.gb1(x1,y1,para)*mx11 + 4*para.gb1(x,y,para)*mx12 + para.gb1(x2,y2,para)*mx13)*normal;
%             e2(i) = e2(i) + w*(para.gb2(x1,y1,para)*mx21 + 4*para.gb2(x,y,para)*mx22 + para.gb2(x2,y2,para)*mx23)*normal;
%         end
    end
end

