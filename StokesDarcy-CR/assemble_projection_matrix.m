function P = assemble_projection_matrix(neighbors,number_of_elements,be,E,dof_u,dof_p,Nx,Ny)

a = 1; b = 0;
a_hat = 0; b_hat = 1;
neighbors_Darcy = neighbors(:,number_of_elements+1:end);

for i = 1:Nx
    for j = 1:Ny
        n1 = (i-1)*2*Ny + 2*j-1; % 当前行列对应的单元
        n2 = n1 + 1; 
        E_new(1,n1) = dof_u + (i-1)*(2*Ny-1) + i*Ny + (j-1)*2 + 1;
        E_new(2,n1) = E_new(1,n1) - (j-1)*2 - (Ny-j) -1;
        E_new(3,n1) = E_new(1,n1) - 1;
        for k = 1:3
            neig = neighbors_Darcy(k,n1);
            if neig == -1
                E_new(k,n1) = -1;
            end
        end
        E_new(1,n2) = E_new(1,n1);
        E_new(2,n2) = E_new(1,n2) + 2*Ny-1 -(j-1);
        E_new(3,n2) = E_new(1,n2) + 1;
        for k = 1:3
            neig = neighbors_Darcy(k,n2);
            if neig == -1   % boundary edge
                E_new(k,n2) = -1;
            end
        end
        
        
        % the right edge is Neumann boundary
%         if i == Nx
%            n2 = (i-1)*2*Ny + 2*j-1 + 1;  
%            if j ~= Ny
%                E_new(1,n2) = dof_u + (i-1)*(2*Ny-1) + i*Ny + (j-1)*2 + 1;
%                E_new(2,n2) = E_new(1,n2) + 2*Ny-1 -(j-1);
%                E_new(3,n2) = E_new(1,n2) + 1;
%            elseif j == Ny
%                E_new(1,n2) = dof_u + (i-1)*(2*Ny-1) + i*Ny + (j-1)*2 + 1;
%                E_new(2,n2) = E_new(1,n2) + 2*Ny-1 -(j-1);
%                E_new(3,n2) = -1;   % 上边界是Dirichlet
%            end
%         end
    end
end
        
nbe = size(be,2);
for k = 1:nbe
    if be(6,k) == -1   % Dirichle 
        edge(k) = be(2,k);  % boudary edge 
    elseif be(6,k) == -2  % Neumann
        edge(k) = 0;
    end
end

number_of_edge_Dirichlet = nnz(edge);  % nnz: compute the number of non-zero element 

Pu = sparse(2*dof_u,2*dof_u- number_of_edge_Dirichlet);
for n = 1:number_of_elements
    for k = 1:3
        neig = neighbors_Darcy(k,n);
        id_ux = E(k,n);   % u1 在当前边的dof编号
        id_uy = E_new(k,n);    % u2
        
        if id_uy ~= -1  % no Dirichlet Boundary
            Pu(id_ux,id_ux) = 1;
            Pu(id_ux+dof_u,id_uy) = 1;  % id_ux+dof_u: 原始uy的编号; id_uy: 处理边界后uy的编号
        else  
            col = find(edge==id_ux);  
            if be(1,col) == 1      % bottom edge
                Pu(id_ux,id_ux) = a;
                Pu(id_ux+dof_u,id_ux) = b;
            elseif be(1,col) == 2  % right edge
                Pu(id_ux,id_ux) = a_hat;
                Pu(id_ux+dof_u,id_ux) = b_hat;
            elseif be(1,col) == 3  % top edge
                Pu(id_ux,id_ux) = a;
                Pu(id_ux+dof_u,id_ux) = b;
            elseif be(1,col) == 4  % left edge
                Pu(id_ux,id_ux) = a_hat;
                Pu(id_ux+dof_u,id_ux) = b_hat;
            end
        end
    end
end
Pp = speye(dof_p);
P = sparse(blkdiag(Pu,Pp));

