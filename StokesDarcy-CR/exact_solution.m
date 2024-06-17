function out = exact_solution(para,Nx_basis,Ny_basis,hx_basis,hy_basis)
global xl xr yb yt
out = zeros((Nx_basis+1)*(Ny_basis+1),1);
for i = 1:Nx_basis+1
    for j = 1:Ny_basis+1
        id = (i-1)*(Ny_basis+1) + j;
        xi = xl + (i-1)*hx_basis;
        yj = yb + (j-1)*hy_basis;
        out(id) = para.u(xi,yj);
    end
end

