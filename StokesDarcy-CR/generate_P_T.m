function [P,T,E,neighbors,neighbors_Stokes,Inter] = generate_P_T(Nx,Ny)
global xl xr yb yt
hx = (xr-xl)/Nx; hy = (yt-yb)/Ny;


%% generate P matrix (保存节点的坐标)
P = zeros(2,(Nx+1)*(Ny+1));
for i = 1:Nx+1
    for j = 1:Ny+1
        id = (i-1)*(Ny+1)+j;  % 全局编号
        xi = xl + (i-1)*hx;   % 当前编号的横坐标
        yj = yb + (j-1)*hy;   % 当前坐标的纵坐标
        P(1,id) = xi;
        P(2,id) = yj;
    end
end
%% generate T matrix (局部单元的全局编号)
for i = 1:Nx
    for j = 1:Ny
        n1 = (i-1)*2*Ny + 2*j-1; %当前行列对应的单元
        n2 = n1 + 1;
        rn1_1 = i;   cn1_1 = j;   % 第n1个单元1号点所在的行和列
        rn1_2 = i+1; cn1_2 = j;   % 第n1个单元2号点所在的行和列
        rn1_3 = i;   cn1_3 = j+1; % 第n1个单元3号点所在的行和列
        jn1_1 = (i-1)*(Ny+1)+j;    % 第n1个单元1号点所对应的全局编号
        jn1_2 = (i)*(Ny+1)+j; 
        jn1_3 = (i-1)*(Ny+1)+j+1; 
        T(1,n1) = jn1_1; T(2,n1) = jn1_2; T(3,n1) = jn1_3;
        rn2_1 = i+1;   rn2_1 = j+1;
        rn2_2 = i;     rn2_2 = j+1;
        rn2_3 = i+1;   rn2_3 = j;
        jn2_1 = (i)*(Ny+1)+j+1;
        jn2_2 = (i-1)*(Ny+1)+j+1;
        jn2_3 = (i)*(Ny+1)+j;
        T(1,n2) = jn2_1; T(2,n2) = jn2_2; T(3,n2) = jn2_3;
    end
end
%% 局部单元边的全局编号
for i = 1:Nx
    for j = 1:Ny
        n1 = (i-1)*2*Ny + 2*j-1;
        n2 = n1+1;
        E(1,n1) = i*Ny + (i-1)*(2*Ny+1) + (j-1)*2 + 1 + 1;
        E(2,n1) = (i-1)*Ny + (i-1)*(2*Ny+1) + j;
        E(3,n1) = i*Ny + (i-1)*(2*Ny+1) + (j-1)*2 + 1;
        E(1,n2) = E(1,n1);
        E(2,n2) = i*Ny + i*(2*Ny+1) + j;
        E(3,n2) = E(1,n2) + 1;
    end
end

%% The following routine finds neighbouring elements in a triangle mesh
np = size(P,2); nt = size(T,2);
% np = (Nx/2+1)*(Ny+1); nt = size(T,2)/2;   % 只生成Omega_f区域的相邻单元
n2e = sparse(np,nt);  % node-to-element adjacency matrix 
                      % n2e(i,j) = 1 means node "i" is in element "j"
for i = 1:nt
    n2e(T(1:3,i),i) = ones(3,1);
end
neighbors = -ones(nt,3);  % -1 means no neighbor
for i=1:nt
    % 1st edge lies between nodes t(2,i) and t(3,i), so search
    % the adjacency matrix for elements sharing these two nodes
    nb=intersect(find(n2e(T(2,i),:)),find(n2e(T(3,i),:)));
    nb=setdiff(nb,i); % remove element "i" from neighbors "nb"
    if isscalar(nb), neighbors(i,1)=nb(1); end
    % 2nd edge
    nb=intersect(find(n2e(T(3,i),:)),find(n2e(T(1,i),:)));
    nb=setdiff(nb,i);
    if isscalar(nb), neighbors(i,2)=nb(1); end
    % 3rd edge
    nb=intersect(find(n2e(T(1,i),:)),find(n2e(T(2,i),:)));
    nb=setdiff(nb,i);
    if isscalar(nb), neighbors(i,3)=nb(1); end
end
neighbors = neighbors';

%% The following routine finds neighbouring elements in a triangle mesh
% np = size(P,2); nt = size(T,2);
np = (Nx/2+1)*(Ny+1); nt = size(T,2)/2;   % 只生成Omega_f区域的相邻单元
n2e = sparse(np,nt);  % node-to-element adjacency matrix 
                      % n2e(i,j) = 1 means node "i" is in element "j"
for i = 1:nt
    n2e(T(1:3,i),i) = ones(3,1);
end
neighbors_Stokes = -ones(nt,3);  % -1 means no neighbor
for i=1:nt
    % 1st edge lies between nodes t(2,i) and t(3,i), so search
    % the adjacency matrix for elements sharing these two nodes
    nb=intersect(find(n2e(T(2,i),:)),find(n2e(T(3,i),:)));
    nb=setdiff(nb,i); % remove element "i" from neighbors "nb"
    if isscalar(nb), neighbors_Stokes(i,1)=nb(1); end
    % 2nd edge
    nb=intersect(find(n2e(T(3,i),:)),find(n2e(T(1,i),:)));
    nb=setdiff(nb,i);
    if isscalar(nb), neighbors_Stokes(i,2)=nb(1); end
    % 3rd edge
    nb=intersect(find(n2e(T(1,i),:)),find(n2e(T(2,i),:)));
    nb=setdiff(nb,i);
    if isscalar(nb), neighbors_Stokes(i,3)=nb(1); end
end
neighbors_Stokes = neighbors_Stokes';

%% 保存界面上的信息
Nx2 = Nx/2;
% Inter.P = zeros(2,Ny+1);  % 界面的顶点坐标
Inter.Tf = zeros(2,Ny);   % T(1,n):第n条界面边所在的局部单元； T(2,n):第n条界面边所在的整体单元

% Inter.P = P(:,Nx*(Ny+1)+1:(Nx+1)*(Ny+1));
for i = 1:Ny
    Inter.Tf(1,i) = (Nx2-1)*Ny*2 + 2*i;  
    Inter.Tf(2,i) = (Nx2-1)*Ny*2 + 2*i;  
end

Inter.Tp = zeros(2,Ny);  % T(1,n):第n条界面边所在的局部单元； T(2,n):第n条界面边所在的整体单元

for i = 1:Ny
    Inter.Tp(1,i) = (2*i-1); % Omega_p 区域单元重新编号
    Inter.Tp(2,i) = Nx2*Ny*2 + (2*i-1);
end

