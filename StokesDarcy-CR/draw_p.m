function draw_p(uh,fun,P,T,number_of_elements,tit)

T = T(:,1:number_of_elements);
dof_p = 1;  % P0 element
figure; 
subplot(1,2,1); hold on;
for ele = 1:number_of_elements
    loc_0 = (1:dof_p)+dof_p*(ele-1);
    uh1 = uh(loc_0);
    trisurf([1 2 3],P(1,T(1:3,ele)),P(2,T(1:3,ele)),full([uh1;uh1;uh1]),'FaceColor','interp','EdgeColor','interp');
end
% t = ones(1,size(T,2)); T = [T;t];
% pdesurf(P,T,uh1);
colormap default; colorbar; view(3);
title('numerical solution of %s', tit)

% 精确解
subplot(1,2,2);hold on;
exact = fun;
for n = 1:number_of_elements
    vertices = P(:,T(1:3,n));
    trisurf([1 2 3],P(1,T(1:3,n)),P(2,T(1:3,n)),...
        eye(3)*[exact(vertices(1,1),vertices(2,1)); exact(vertices(1,2),vertices(2,2));exact(vertices(1,3),vertices(2,3))],...
        'FaceColor','interp','EdgeColor','interp');
end
title('exact solution of %s', tit); view(3);
colorbar;


% figure;hold on;
% plot(uh,'b');
% for n = 1:number_of_elements
%     vertices = P(:,T(:,n));
%     x = (vertices(1,1) + vertices(1,2) + vertices(1,3))/3;
%     y = (vertices(2,1) + vertices(2,2) + vertices(2,3))/3;
%     p(n,1) = fun(x,y);
% end
% plot(p,'r'); 
% legend('numerical','exact');

