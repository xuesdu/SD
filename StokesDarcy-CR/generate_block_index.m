function [blk_Stokes, blk_Darcy, Blk_Stokes, Blk_Darcy] = generate_block_index(T, E, neighbors, Nx, Ny, dof_us)

number_of_nodes = (Nx+1)*(Ny+1);
number_of_elements_Stokes = size(T,2)/2;
T_Stokes = T(:,1:number_of_elements_Stokes);

for k = 1:number_of_nodes
    [row,col] = find(T_Stokes == k);   
    l = length(col);
    Row(1:l,k) = row;
    Col(1:l,k) = col;
end
%%% ----------------- us1 -------------------------------------
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
    As1(1:l,j) = temp;
end
blk_us1 = sort(As1);
%%%% --------------------- us2 --------------------------------------
for i = 1:size(As1,1)
    for j = 1:size(As1,2)
        if As1(i,j) == 0
            As2(i,j) = 0;
        else
            As2(i,j) = As1(i,j) + dof_us;
        end
    end
end
blk_us2 = sort(As2);
blk_Stokes = [blk_us1;blk_us2];

%%%% ---------------------- Darcy domain ----------------------------------
neighbors_Darcy = neighbors(:,number_of_elements_Stokes+1:end);
for i = 1:Nx
    for j = 1:Ny
        n1 = (i-1)*2*Ny + 2*j-1; % 当前行列对应的单元
        n2 = n1 + 1; 
        %%% ud1 
        E_ud1(1,n1) = (i-1)*(2*Ny+1) + (i-1)*Ny + 2*j;
        if i == 1
            E_ud1(2,n1) = 0;
        else
            E_ud1(2,n1) = (i-2)*Ny + (i-1)*(2*Ny+1) + j;
        end
        E_ud1(3,n1) = E_ud1(1,n1) - 1;
        
        E_ud1(1,n2) = E_ud1(1,n1);
        E_ud1(2,n2) = (i-1)*Ny + i*(2*Ny+1) + j;
        E_ud1(3,n2) = E_ud1(1,n2) + 1;
        
        %%% ud2
        E_ud2(1,n1) = dof_us + (i-1)*(2*Ny-1) + i*Ny + (j-1)*2 + 1 - Ny;
        E_ud2(2,n1) = E_ud2(1,n1) - (j-1)*2 - (Ny-j) -1;
        E_ud2(3,n1) = E_ud2(1,n1) - 1;
        for k = 1:3
            neig = neighbors_Darcy(k,n1);
            if neig == -1
                E_ud2(k,n1) = 0;
            end
        end
        E_ud2(1,n2) = E_ud2(1,n1);
        E_ud2(2,n2) = E_ud2(1,n2) + 2*Ny-1 -(j-1);
        E_ud2(3,n2) = E_ud2(1,n2) + 1;
        for k = 1:3
            neig = neighbors_Darcy(k,n2);
            if neig == -1   % boundary edge
                E_ud2(k,n2) = 0;
            end
        end
    end
end

%%%% ----------------------- ud1 ---------------------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            AA(2*i,j) = 0;
            AA(2*i-1,j) = 0;
        else
            if Row(i,j) == 1
                AA(2*i,j) = E_ud1(2,Col(i,j));
                AA(2*i-1,j) = E_ud1(3,Col(i,j));
            elseif Row(i,j) == 2
                AA(2*i,j) = E_ud1(3,Col(i,j));
                AA(2*i-1,j) = E_ud1(1,Col(i,j));
            elseif Row(i,j) == 3
                AA(2*i,j) = E_ud1(1,Col(i,j));
                AA(2*i-1,j) = E_ud1(2,Col(i,j));
            end
        end
    end
end

for j = 1:size(AA,2)
    temp = unique(AA(:,j));
    l = length(temp);
    Ad1(1:l,j) = temp;
end
% add the dofs of Stokes velocity 
for i = 1:size(Ad1,1)
    for j = 1:size(Ad1,2)
        if Ad1(i,j) == 0
            Ad1(i,j) = 0;
        else
            Ad1(i,j) = Ad1(i,j) + 2*dof_us;
        end
    end
end
blk_ud1 = sort(Ad1);
%%%% ----------------------- ud2 ---------------------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            AA(2*i,j) = 0;
            AA(2*i-1,j) = 0;
        else
            if Row(i,j) == 1
                AA(2*i,j) = E_ud2(2,Col(i,j));
                AA(2*i-1,j) = E_ud2(3,Col(i,j));
            elseif Row(i,j) == 2
                AA(2*i,j) = E_ud2(3,Col(i,j));
                AA(2*i-1,j) = E_ud2(1,Col(i,j));
            elseif Row(i,j) == 3
                AA(2*i,j) = E_ud2(1,Col(i,j));
                AA(2*i-1,j) = E_ud2(2,Col(i,j));
            end
        end
    end
end
for j = 1:size(AA,2)
    temp = unique(AA(:,j));
    l = length(temp);
    Ad2(1:l,j) = temp;
end
% add the dofs of Stokes velocity 
for i = 1:size(Ad2,1)
    for j = 1:size(Ad2,2)
        if Ad2(i,j) == 0
            Ad2(i,j) = 0;
        else
            Ad2(i,j) = Ad2(i,j) + 2*dof_us;
        end
    end
end
blk_ud2 = sort(Ad2);
blk_Darcy = [blk_ud1;blk_ud2];

%%%% 顶点所在单元涉及的所有自由度
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            BB(3*i-2:3*i,j) = 0;
        else
            BB(3*i-2:3*i,j) = E(:,Col(i,j));
        end
    end
end
%%%%% ------------ us1 ----------------
for j = 1:size(BB,2)
    temp = unique(BB(:,j));
    l = length(temp);
    Bs1(1:l,j) = temp;
end
Blk_us1 = sort(Bs1);
%%%%% ------------ us2 ------------------
for i = 1:size(Bs1,1)
    for j = 1:size(Bs1,2)
        if Bs1(i,j) == 0
            Bs2(i,j) = 0;
        else
            Bs2(i,j) = Bs1(i,j) + dof_us;
        end
    end
end
Blk_us2 = sort(Bs2);
Blk_Stokes = [Blk_us1;Blk_us2];
%%%% ------------ ud1 -------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            BB(3*i-2:3*i,j) = 0;
        else
            BB(3*i-2:3*i,j) = E_ud1(:,Col(i,j));
        end
    end
end
for j = 1:size(BB,2)
    temp = unique(BB(:,j));
    l = length(temp);
    Bd1(1:l,j) = temp;
end
% add the dofs of Stokes velocity 
for i = 1:size(Bd1,1)
    for j = 1:size(Bd1,2)
        if Bd1(i,j) == 0
            Bd1(i,j) = 0;
        else
            Bd1(i,j) = Bd1(i,j) + 2*dof_us;
        end
    end
end
Blk_ud1 = sort(Bd1);
%%%% ------------ ud2 -------------------
for i = 1:size(Col,1)
    for j = 1:size(Col,2)
        if Row(i,j) == 0
            BB(3*i-2:3*i,j) = 0;
        else
            BB(3*i-2:3*i,j) = E_ud2(:,Col(i,j));
        end
    end
end
for j = 1:size(BB,2)
    temp = unique(BB(:,j));
    l = length(temp);
    Bd2(1:l,j) = temp;
end
% add the dofs of Stokes velocity 
for i = 1:size(Bd2,1)
    for j = 1:size(Bd2,2)
        if Bd2(i,j) == 0
            Bd2(i,j) = 0;
        else
            Bd2(i,j) = Bd2(i,j) + 2*dof_us;
        end
    end
end
Blk_ud2 = sort(Bd2);

Blk_Darcy = [Blk_ud1;Blk_ud2];

% make them sparse -- XH
blk_Stokes = sparse(blk_Stokes);
blk_Darcy = sparse(blk_Darcy);
Blk_Stokes = sparse(Blk_Stokes);
Blk_Darcy = sparse(Blk_Darcy);

end