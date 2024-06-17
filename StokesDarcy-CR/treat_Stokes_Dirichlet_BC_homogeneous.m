function [A,b] = treat_Stokes_Dirichlet_BC_homogeneous(A,b,be_s,P,dof_us)

nbe_s = size(be_s,2);

for k = 1:nbe_s
    if be_s(6,k) == -1 % Dirichlet        
        edge_index = be_s(2,k);
        % u1
        A(edge_index, :) = 0; 
        A(:, edge_index) = 0;
        A(edge_index, edge_index) = 1; 
        b(edge_index,1) = 0;
        % u2
        A(edge_index+dof_us, :) = 0; 
        A(:, edge_index + dof_us) = 0;
        A(edge_index+dof_us, edge_index+dof_us) = 1; 
        b(edge_index+dof_us,1) = 0;
    end
end


