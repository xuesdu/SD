function draw(u1,u2,para,hx_basis_u,hy_basis_u)

global xl xr yb yt xbar
Nx_u = (xbar-xl)/hx_basis_u; Ny_u = (yt-yb)/hy_basis_u;
for i = 1:Nx_u+1
    for j = 1:Ny_u+1
        xi = xl + (i-1)*hx_basis_u;
        yj = yb + (j-1)*hy_basis_u;
        id = (i-1)*(Ny_u+1)+j;  % 全局编号
        u1_exact(id) = para.us1(xi,yj);
        u2_exact(id) = para.us2(xi,yj);
    end
end
[Xu,Yu] = meshgrid(xl:hx_basis_u:xbar,yb:hy_basis_u:yt);
uh1 = reshape(u1,Nx_u+1,Ny_u+1); U1 = reshape(u1_exact,Nx_u+1,Ny_u+1); e1 = uh1-U1;
uh2 = reshape(u2,Nx_u+1,Ny_u+1); U2 = reshape(u2_exact,Nx_u+1,Ny_u+1); e2 = uh2-U2;
figure(1); subplot(1,3,1); surf(Xu,Yu,uh1); subplot(1,3,2); surf(Xu,Yu,U1); subplot(1,3,3); surf(Xu,Yu,e1); 
figure(2); subplot(1,3,1); surf(Xu,Yu,uh2); subplot(1,3,2); surf(Xu,Yu,U2); subplot(1,3,3); surf(Xu,Yu,e2); 

