function draw_p(uh,pp,para,Nx,Ny)
global number_of_elements P T 
P_s = P(:,1:(Nx+1)*(Ny+1));
T_s = T(:,1:number_of_elements/2);
dof_p = 1;  % P0 element
figure; subplot(1,2,1);hold on;
for n = 1:number_of_elements/2
    loc_0 = (1:dof_p)+dof_p*(n-1);
    uh1(n) = uh(loc_0);
end
t = ones(1,size(T_s,2)); T_s = [T_s;t];
pdesurf(P_s,T_s,uh1);colormap default; colorbar; view(3);
title('numerical solution of pf')
subplot(1,2,2); hold on;

% 精确解
exact = para.ps;
for n = 1:number_of_elements/2
    vertices = P_s(:,T_s(1:3,n));
    trisurf([1 2 3],P_s(1,T_s(1:3,n)),P_s(2,T_s(1:3,n)),...
        eye(3)*[exact(vertices(1,1),vertices(2,1)); exact(vertices(1,2),vertices(2,2));exact(vertices(1,3),vertices(2,3))],...
        'FaceColor','interp','EdgeColor','interp');
end
colorbar;
title('exact solution of pf');view(3);


% numerical solution of pd
P_d = P(:,Nx*(Ny+1)+1:end);
T_d = T(:,number_of_elements/2+1:end);
dof_p = 1;  % P0 element
figure; subplot(1,2,1); hold on;
for n = 1:number_of_elements/2
    loc_0 = (1:dof_p)+dof_p*(n-1);
    uh1(n) = pp(loc_0);
end
t = ones(1,size(T_d,2)); T_d = [T_d;t];
pdesurf(P_d,T_s,uh1);colormap default; colorbar; view(3);
title('numerical solution of pp')
subplot(1,2,2); hold on;
% exact solution of pd
exact = para.pd;
for n = 1:number_of_elements/2
    vertices = P(:,T_d(1:3,n));
    trisurf([1 2 3],P(1,T_d(1:3,n)),P(2,T_d(1:3,n)),...
        eye(3)*[exact(vertices(1,1),vertices(2,1)); exact(vertices(1,2),vertices(2,2));exact(vertices(1,3),vertices(2,3))],...
        'FaceColor','interp','EdgeColor','interp');
end
colorbar;
title('exact solution of pp'); view(3);

