function [bn,be] = generate_boundary_nodes_edges(Nx_basis,Ny_basis,Nx,Ny)

nbn = 2*(Nx_basis+Ny_basis);
bn = zeros(2,nbn);
bn(1,:) = -1;  % Dirichlet boundary condition
% 下边界
for i = 1:Nx_basis+1
    k = i;      % 第 k 个边界点
    id = (i-1)*(Ny_basis+1) + 1;  % 全局编号
    bn(2,k) = id;
end
% 右边界
for j = 1:Ny_basis
    k = Nx_basis+j;
    id = Nx_basis*(Ny_basis+1) + j;
    bn(2,k) = id;
end
% 上边界
for i = 1:Nx_basis
    k = Nx_basis+Ny_basis+i;
    id = (Nx_basis+1)*(Ny_basis+1) - (i-1)*(Ny_basis+1);
    bn(2,k) = id;
end
% 左边界
for j = 1:Ny_basis
    k = 2*Nx_basis+Ny_basis+j;
    id = Ny_basis+1-(j-1);
    bn(2,k) = id;
end

nbe = 2*(Nx+Ny);
be = zeros(5,nbe);
% be(1,k): on which edge of Omega 1,2,3,4: 1-bottom 2-right 3-top 4-left
% be(2,k): the index of the edge
% be(3,k): the index of the element which contains the kth boundary edge
% be(3,k): the global node index of the first end node
% be(4,k): the global node index of the second end node

% 下边界
for i = 1:Nx
    k = i;  % the kth boundary edge
    be(1,k) = 1;  % bottom edge
    edge = i*Ny + (i-1)*(2*Ny+1) + 1; % edge index
    element = (i-1)*Ny*2 + 1;  % element index
    nl = (i-1)*(Ny+1) + 1;
    nr = i*(Ny+1) + 1;
    be(2,k) = edge;
    be(3,k) = element;
    be(4,k) = nl;
    be(5,k) = nr;
end
% 右边界
for j = 1:Ny
    k = Nx + j;
    be(1,k) = 2; 
    edge = Nx*Ny + Nx*(2*Ny+1) + j;
    element = (Nx-1)*Ny*2 + j*2;
    nl = Nx*(Ny+1) + j;
    nr = Nx*(Ny+1) + j+1;
    be(2,k) = edge;
    be(3,k) = element;
    be(4,k) = nl;
    be(5,k) = nr;
end
% 上边界
for i = 1:Nx
    k = Nx + Ny + i;
    be(1,k) = 3;
    edge = Nx*Ny + Nx*(2*Ny+1) - (i-1)*(2*Ny+1) - (i-1)*Ny;
    element = 2*Nx*Ny - (i-1)*2*Ny;
    nl = (Nx+1)*(Ny+1) - (i-1)*(Ny+1);
    nr = (Nx+1)*(Ny+1) - (i)*(Ny+1);
    be(2,k) = edge;
    be(3,k) = element;
    be(4,k) = nl;
    be(5,k) = nr;
end
% 左边界
for j = 1:Ny
    k = 2*Nx + Ny + j;
    be(1,k) = 4;
    edge = Ny - (j-1);
    element = 2*Ny-(2*j-1);
    nl = Ny+1-(j-1);
    nr = Ny+1-j;
    be(2,k) = edge;
    be(3,k) = element;
    be(4,k) = nl;
    be(5,k) = nr;
end
