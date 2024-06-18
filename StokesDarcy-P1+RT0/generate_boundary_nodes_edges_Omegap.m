function [bn,be] = generate_boundary_nodes_edges_Omegap(Nx_basis,Ny_basis,Nx,Ny,BC)

nbn = 2*(Nx_basis+Ny_basis);
bn = zeros(2,nbn);
% 下边界
for i = 1:Nx_basis
    k = i;      % 第 k 个边界点
    id = (i-1)*(Ny_basis+1) + 1;  % 全局编号
    bn(2,k) = id;
    switch BC
        case 'I' 
            bn(1,k) = -2;
        case 'II'
            bn(1,k) = -1;
        case 'III'
            bn(1,k) = -1;
        case 'IV'
            bn(1,k) = -2;
    end
end
% 右边界
for j = 1:Ny_basis
    k = Nx_basis+j;
    id = Nx_basis*(Ny_basis+1) + j;
    bn(2,k) = id;
    switch BC
        case 'I' 
            bn(1,k) = -2;
        case 'II'
            bn(1,k) = -1;
        case 'III'
            bn(1,k) = -1;
        case 'IV'
            bn(1,k) = -2;
    end
end
% 上边界
for i = 1:Nx_basis+1
    k = Nx_basis+Ny_basis+i;
    id = (Nx_basis+1)*(Ny_basis+1) - (i-1)*(Ny_basis+1);
    bn(2,k) = id;
    switch BC
        case 'I' 
            bn(1,k) = -2;
        case 'II'
            bn(1,k) = -1;
        case 'III'
            bn(1,k) = -1;
        case 'IV'
            bn(1,k) = -2;
    end
end
% 左边界
% % for j = 1:Ny_basis
% %     k = 2*Nx_basis+Ny_basis+j;
% %     id = Ny_basis+1-(j-1);
% %     bn(2,k) = id;
% % end

nbe = 2*(Nx+Ny);
be = zeros(6,nbe);
% be(1,k): on which edge of Omega 1,2,3,4: 1-bottom 2-right 3-top 4-left
% be(2,k): the index of the edge
% be(3,k): the index of the element which contains the kth boundary edge
% be(4,k): the global node index of the first end node
% be(5,k): the global node index of the second end node
% be(6,k): the type of boudary condition

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
    switch BC
        case 'I' 
            be(6,k) = -2;
        case 'II'
            be(6,k) = -1;
        case 'III'
            be(6,k) = -1;
        case 'IV'
            be(6,k) = -2;
    end
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
    switch BC
        case 'I' 
            be(6,k) = -2;
        case 'II'
            be(6,k) = -1;
        case 'III'
            be(6,k) = -1;
        case 'IV'
            be(6,k) = -2;
    end
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
    switch BC
        case 'I' 
            be(6,k) = -2;
        case 'II'
            be(6,k) = -1;
        case 'III'
            be(6,k) = -1;
        case 'IV'
            be(6,k) = -2;
    end
end
% 左边界
% % for j = 1:Ny
% %     k = 2*Nx + Ny + j;
% %     be(1,k) = 4;
% %     edge = Ny - (j-1);
% %     element = 2*Ny-(2*j-1);
% %     nl = Ny+1-(j-1);
% %     nr = Ny+1-j;
% %     be(2,k) = edge;
% %     be(3,k) = element;
% %     be(4,k) = nl;
% %     be(5,k) = nr;
% % end
