function ug = generate_nonhomogeneous_BC(Omega,para,bn,dof_u)
global P
if Omega == 's'  % P1
    nbn_s = size(bn,2);
    ug = sparse(2*dof_u,1);
    for k = 1:nbn_s
        if bn(1,k) == -1 
            id = bn(2,k);
            x = P(1,id); y = P(2,id);  
            u1_exact = para.us1(x,y);
            u2_exact = para.us2(x,y);
            ug(id) = u1_exact;
            ug(id+dof_u) = u2_exact;
        end
    end  
end
