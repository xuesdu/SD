function error = Error_L2_p(Omega,solution,fun,Gpn,hx,hy)
global number_of_elements P T
ne = number_of_elements/2;
Gauss_weights = 2*[-27/48;25/48;25/48;25/48];
Gauss_nodes = [1/3 3/5 1/5 1/5;1/3 1/5 1/5 3/5];
error = 0;
for n = 1:ne
    ph_local = solution(n);  % 当前单元上的数值解
    if Omega == 1
        vertices = P(:,T(:,n));
    elseif Omega == 2
        vertices = P(:,T(:,n+number_of_elements/2));
    end
    x1 = vertices(1,1); y1 = vertices(2,1);
    x2 = vertices(1,2); y2 = vertices(2,2);
    x3 = vertices(1,3); y3 = vertices(2,3);
    xn = [x1 x2 x3]; yn = [y1 y2 y3];
    Jacobi = [xn(2)-xn(1) xn(3)-xn(1);yn(2)-yn(1) yn(3)-yn(1)];
    
    p = Jacobi*Gauss_nodes+[xn(1) 0;yn(1) 0]*ones(2,Gpn);
    error_temp = 0;
    for i = 1:Gpn
        exact = fun(p(1,i),p(2,i));
        error_temp = error_temp + Gauss_weights(i)*((exact-ph_local)^2);
    end
    error = error + 0.5*hx*hy*error_temp;
end

error = sqrt(error);