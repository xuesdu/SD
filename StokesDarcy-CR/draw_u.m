function draw_u(uh,fun,P,T,E,number_of_elements,tit)
dof_u = 3; % CR元每个单元上的自由度
figure;
subplot(1,2,1); hold on;
for n = 1:number_of_elements
    uh1 = uh(E(:,n));
    trisurf([1 2 3],P(1,T(1:3,n)),P(2,T(1:3,n)),full(uh1),'FaceColor','interp','EdgeColor','interp');
end
title(sprintf('numerical solution of %s', tit));colorbar; view(3);

subplot(1,2,2); hold on;
for i = 1:number_of_elements
    vtx = P(:,T(1:3,i));
    trisurf([1 2 3],P(1,T(1:3,i)),P(2,T(1:3,i)),...
        eye(3)*[fun(vtx(1,1),vtx(2,1));fun(vtx(1,2),vtx(2,2));fun(vtx(1,3),vtx(2,3))],...
        'FaceColor','interp','EdgeColor','interp');
end
title(sprintf('exact solution of %s', tit));colorbar;view(3);

