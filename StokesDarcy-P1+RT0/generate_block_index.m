function [blk_Stokes, blk_Darcy, Blk_Stokes, Blk_Darcy] = generate_block_index(Nx, Ny, dof_ul, dof_uR)

global T E
number_of_nodes_Stokes = (Nx+1)*(Ny+1);
number_of_elements_Stokes = size(T,2)/2;
T_Stokes = T(:,1:number_of_elements_Stokes);

for k = 1:number_of_nodes_Stokes
    [row,col] = find(T_Stokes == k);   
    l = length(col);
    Row(1:l,k) = row;
    Col(1:l,k) = col;
end
%%% -------------------- us1 (P1) -------------------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            AA(2*i,j) = 0;
            AA(2*i-1,j) = 0;
        else
            if Row(i,j) == 1
                AA(2*i,j) = T(2,Col(i,j));
                AA(2*i-1,j) = T(3,Col(i,j));
            elseif Row(i,j) == 2
                AA(2*i,j) = T(3,Col(i,j));
                AA(2*i-1,j) = T(1,Col(i,j));
            elseif Row(i,j) == 3
                AA(2*i,j) = T(1,Col(i,j));
                AA(2*i-1,j) = T(2,Col(i,j));
            end
        end
    end
end

for j = 1:size(AA,2)
    temp = unique(AA(:,j));
    l = length(temp);
    As1(1:l,j) = temp;
end
blk_us1 = sort(As1);
%%%% --------------------- us2 (P1) --------------------------------------
for i = 1:size(As1,1)
    for j = 1:size(As1,2)
        if As1(i,j) == 0
            As2(i,j) = 0;
        else
            As2(i,j) = As1(i,j) + dof_ul;
        end
    end
end
blk_us2 = sort(As2);
%%%% -------------------- usR (RT0) ------------------------------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            AA(2*i,j) = 0;
            AA(2*i-1,j) = 0;
        else
            if Row(i,j) == 1
                AA(2*i,j) = E(2,Col(i,j));
                AA(2*i-1,j) = E(3,Col(i,j));
            elseif Row(i,j) == 2
                AA(2*i,j) = E(3,Col(i,j));
                AA(2*i-1,j) = E(1,Col(i,j));
            elseif Row(i,j) == 3
                AA(2*i,j) = E(1,Col(i,j));
                AA(2*i-1,j) = E(2,Col(i,j));
            end
        end
    end
end
for j = 1:size(AA,2)
    temp = unique(AA(:,j));
    l = length(temp);
    AsR(1:l,j) = temp;
end
% add the dofs of usl
for i = 1:size(AsR,1)
    for j = 1:size(AsR,2)
        if AsR(i,j) == 0
            AsR(i,j) = 0;
        else
            AsR(i,j) = AsR(i,j) + 2*dof_ul;
        end
    end
end
blk_usR = sort(AsR);

blk_Stokes = [blk_us1;blk_us2;blk_usR];

%%%% Darcy domain
for i = 1:Nx
    for j = 1:Ny
        n1 = (i-1)*2*Ny + 2*j-1; % 当前行列对应的单元
        n2 = n1 + 1; 
        %%% ud1 
        E_ud(1,n1) = (i-1)*(2*Ny+1) + (i-1)*Ny + 2*j;
        if i == 1
            E_ud(2,n1) = 0;
        else
            E_ud(2,n1) = (i-2)*Ny + (i-1)*(2*Ny+1) + j;
        end
        E_ud(3,n1) = E_ud(1,n1) - 1;
        
        E_ud(1,n2) = E_ud(1,n1);
        E_ud(2,n2) = (i-1)*Ny + i*(2*Ny+1) + j;
        E_ud(3,n2) = E_ud(1,n2) + 1;
        
    end
end
%%%% --------------- ud (RT0) ------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            AA(2*i,j) = 0;
            AA(2*i-1,j) = 0;
        else
            if Row(i,j) == 1
                AA(2*i,j) = E_ud(2,Col(i,j));
                AA(2*i-1,j) = E_ud(3,Col(i,j));
            elseif Row(i,j) == 2
                AA(2*i,j) = E_ud(3,Col(i,j));
                AA(2*i-1,j) = E_ud(1,Col(i,j));
            elseif Row(i,j) == 3
                AA(2*i,j) = E_ud(1,Col(i,j));
                AA(2*i-1,j) = E_ud(2,Col(i,j));
            end
        end
    end
end
for j = 1:size(AA,2)
    temp = unique(AA(:,j));
    l = length(temp);
    AdR(1:l,j) = temp;
end
% add the dofs of usl
for i = 1:size(AdR,1)
    for j = 1:size(AdR,2)
        if AdR(i,j) == 0
            AdR(i,j) = 0;
        else
            AdR(i,j) = AdR(i,j) + 2*dof_ul + dof_uR;
        end
    end
end
blk_Darcy = sort(AdR);

% ---------------------------------------------------------------------------------------
%%%% 顶点所在单元涉及的所有自由度

%%%%% ------------ us1 (P1) ----------------
Blk_us1 = blk_us1;
%%%%% ------------ us2 (P2) ------------------
Blk_us2 = blk_us2;
%%%% ------------ usR(RT0) -------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            BB(3*i-2:3*i,j) = 0;
        else
            BB(3*i-2:3*i,j) = E(:,Col(i,j));
        end
    end
end
for j = 1:size(BB,2)
    temp = unique(BB(:,j));
    l = length(temp);
    BsR(1:l,j) = temp;
end
% add the dofs of usl
for i = 1:size(BsR,1)
    for j = 1:size(BsR,2)
        if BsR(i,j) == 0
            BsR(i,j) = 0;
        else
            BsR(i,j) = BsR(i,j) + 2*dof_ul;
        end
    end
end
Blk_usR = sort(BsR);
Blk_Stokes = [Blk_us1; Blk_us2; Blk_usR];
%%%% ------------ ud1 -------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            BB(3*i-2:3*i,j) = 0;
        else
            BB(3*i-2:3*i,j) = E_ud(:,Col(i,j));
        end
    end
end
for j = 1:size(BB,2)
    temp = unique(BB(:,j));
    l = length(temp);
    BdR(1:l,j) = temp;
end
% add the dofs of Stokes velocity 
for i = 1:size(BdR,1)
    for j = 1:size(BdR,2)
        if BdR(i,j) == 0
            BdR(i,j) = 0;
        else
            BdR(i,j) = BdR(i,j) + 2*dof_ul + dof_uR;
        end
    end
end
Blk_Darcy = sort(BdR);


end
