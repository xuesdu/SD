function [P,T,E,Inter] = generate_P_T(Nx2,Ny,basis_type)
global xl xr yb yt
hx = (xr-xl)/Nx2; hy = (yt-yb)/Ny;

if basis_type == 201
    % generate P matrix (保存节点的坐标)
    P = zeros(2,(Nx2+1)*(Ny+1));
    for i = 1:Nx2+1
        for j = 1:Ny+1
            id = (i-1)*(Ny+1)+j;  % 全局编号
            xi = xl + (i-1)*hx;   % 当前编号的横坐标
            yj = yb + (j-1)*hy;   % 当前坐标的纵坐标
            P(1,id) = xi;
            P(2,id) = yj;
        end
    end
    
    % generate T matrix (局部单元的全局编号)
    for i = 1:Nx2
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
            rn2_1 = i;     rn2_1 = j+1;
            rn2_2 = i+1;   rn2_2 = j;
            rn2_3 = i+1;   rn2_3 = j+1;
            jn2_1 = (i)*(Ny+1)+j+1;
            jn2_2 = (i-1)*(Ny+1)+j+1;
            jn2_3 = (i)*(Ny+1)+j;
            T(1,n2) = jn2_1; T(2,n2) = jn2_2; T(3,n2) = jn2_3;
        end
    end
    % E: 局部单元边的全局编号
    for i = 1:Nx2
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
    
    
elseif basis_type == 202
    hx2 = hx/2; hy2 = hy/2;
    Nx2 = Nx2*2; Ny2 = Ny*2;
    % 生成 P 矩阵
    for i = 1:Nx2+1
        for j = 1:Ny2+1
            xi = xl + (i-1)*hx2;
            yj = yb + (j-1)*hy2;
            id = (i-1)*(Ny2+1) + j;
            P(1,id) = xi;
            P(2,id) = yj;
        end
    end
    % generate T matrix (局部单元的全局编号)
    % 所有的编号存到Q矩阵中
    for i = 1:Nx2+1
        for j = 1:Ny2+1
            Q(i,j) = (i-1)*(Ny2+1) + j;
        end
    end
    % 找出所在的行列
    for n = 1:Nx2*Ny
        if mod(n,Ny) == 0
            row = Ny;
            column = n/Ny;
        else
            row = mod(n,Ny);
            column = fix(n/Ny)+1; % fix: 向下取整
        end
        
        n1 = 2*n-1;  % 所对应的单元
        n2 = 2*n;
        T(1,n1) = Q(2*column-1,2*row-1);
        T(2,n1) = Q(2*column+1,2*row-1);
        T(3,n1) = Q(2*column-1,2*row+1);
        T(4,n1) = Q(2*column,2*row-1);
        T(5,n1) = Q(2*column,2*row);
        T(6,n1) = Q(2*column-1,2*row);
        
        T(1,n2) = Q(2*column-1,2*row+1);
        T(2,n2) = Q(2*column+1,2*row-1);
        T(3,n2) = Q(2*column+1,2*row+1);
        T(4,n2) = Q(2*column,2*row);
        T(5,n2) = Q(2*column+1,2*row);
        T(6,n2) = Q(2*column,2*row+1);
    end
        
end

% 保存界面上的信息
Nx = Nx2/2;
% Inter.Pf = zeros(2,Ny+1);  
Inter.Tf = zeros(2,Ny);  % T(1,n):第n条界面边所在的单元；T(2,n):第n条界面边所在的下端点； T(3,n):第n条界面边所在的上端点
% Inter.Ef = zeros(1,Ny);

% Inter.Pf = P(:,Nx*(Ny+1)+1:(Nx+1)*(Ny+1));
for i = 1:Ny
    Inter.Tf(1,i) = (Nx-1)*Ny*2 + 2*i;   % 界面边在Stokes区域局部单元编号
    Inter.Tf(2,i) = (Nx-1)*Ny*2 + 2*i;   % 界面边global编号
%     Inter.Tf(3,i) = Nx*(Ny+1) + i + 1;
%     Inter.Ef(1,i) = Nx*Ny + Nx*(2*Ny+1) + i;
end

% Inter.Pp = zeros(2,Ny+1);
Inter.Tp = zeros(2,Ny);  % T(1,n):第n条界面边所在单元的局部编号（Darcy区域从头编号）；T(2,n):第n条界面边所在单元的整体编号
% Inter.Ep = zeros(1,Ny);

% Inter.Pp = P(:,Nx*(Ny+1)+1:(Nx+1)*(Ny+1));
for i = 1:Ny
    Inter.Tp(1,i) = (2*i-1); % Omega_p 区域单元重新编号
    Inter.Tp(2,i) = Nx*Ny*2 + (2*i-1);
%     Inter.Tp(3,i) = i + 1;
end

    
