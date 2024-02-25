clear; clc;
xl = 0; xr = 2; yb = 0; yt = 1; xbar = 1;
Nx = 4; Ny = 2;
hx = (xr-xl)/Nx; hy = (yt-yb)/Ny;
ne = 2*Nx*Ny; % number of elements

% generate P matrix (保存节点的坐标)
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

% generate T matrix (局部单元的全局编号)
T = zeros(3,ne);
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
        rn2_1 = i+1;     rn2_1 = j+1;
        rn2_2 = i;       rn2_2 = j+1;
        rn2_3 = i+1;     rn2_3 = j;
        jn2_1 = (i)*(Ny+1)+j+1;
        jn2_2 = (i-1)*(Ny+1)+j+1;
        jn2_3 = (i)*(Ny+1)+j;
        T(1,n2) = jn2_1; T(2,n2) = jn2_2; T(3,n2) = jn2_3;
    end
end

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

for n = 1:ne
    x1 = P(1,T(1,n)); % 当前单元1号节点的横坐标
    y1 = P(2,T(1,n));
    x2 = P(1,T(2,n)); % 当前单元2号节点的横坐标
    y2 = P(2,T(2,n));
    x3 = P(1,T(3,n));
    y3 = P(2,T(3,n)); % 当前单元3号节点的横坐标
    plot([x1 x2 x3 x1],[y1 y2 y3 y1],'b')
    hold on
    text(x1+hx/100,y1+hy/100,num2str(T(1,n)),'color','r');  % 顶点编号
    text(x2+hx/100,y2+hy/100,num2str(T(2,n)),'color','r');
    text(x3+hx/100,y3+hy/100,num2str(T(3,n)),'color','r');
    if mod(n,2) == 1
        text(x2-hx/2+hx/100,y2+hy/2+hy/100,num2str(E(1,n)),'color','k');  % 边的编号
        text(x3+hx/100,y3-hy/2,num2str(E(2,n)),'color','k');
        text(x1+hx/2,y1+hy/100,num2str(E(3,n)),'color','k');
    else
        text(x2+hx/2+hx/100,y2-hy/2+hy/100,num2str(E(1,n)),'color','k');
        text(x3+hx/100,y3+hy/2,num2str(E(2,n)),'color','k');
        text(x1-hx/2,y1+hy/100,num2str(E(3,n)),'color','k');
    end
    text((x1+x2+x3)/3,(y1+y2+y3)/3,num2str(n),'color','m');  % 单元的编号
end
axis off
        
    