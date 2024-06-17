function [A,bs] = treat_Stokes_Dirichlet_BC(A,bs,para,be_s,P,dof_us)

nbe_s = size(be_s,2);

for k = 1:nbe_s
    if be_s(1,k) ~= 0        
        edge_index = be_s(2,k);
        x1 = P(1,be_s(4,k)); y1 = P(2,be_s(4,k));
        x2 = P(1,be_s(5,k)); y2 = P(2,be_s(5,k));
        x = (x1+x2)/2; y = (y1+y2)/2;
        % u1
        A(edge_index,:) = 0; 
        A(edge_index,edge_index) = 1; 
        bs(edge_index) =  para.us1(x,y);
        % u2
        A(edge_index+dof_us,:) = 0; 
        A(edge_index+dof_us,edge_index+dof_us) = 1; 
        bs(edge_index+dof_us) = para.us2(x,y);
    end
end


