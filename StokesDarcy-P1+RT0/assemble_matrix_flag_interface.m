function PI = assemble_matrix_flag_interface(Ns, Nd, Ny, dof_ul, dof_uR)
% N = Ns+Nd: total dofs; 
% Ny: number of interface edges;

global Inter T E

N = Ns+Nd;
I = eye(Ns+Nd);

dof_vec = [];

for n = 1:Ny
    ele_s = Inter.Tf(1,n);
    ele_d = Inter.Tp(1,n);
    
    % total dofs near the interface
    dof_inter = [T(:,ele_s); T(:,ele_s)+dof_ul; E(:,ele_s)+2*dof_ul; 2*dof_ul+dof_uR+ele_s;...
        Ns+E(1,ele_d)-Ny; Ns+E(3,ele_d)-Ny; Ns+dof_uR-Ny+ele_d];
    % [us1,ud,pd]
    % dof_inter = [T(:,ele_s); Ns+E(1,ele_d)-Ny; Ns+E(3,ele_d)-Ny; Ns+dof_uR-Ny+ele_d];

    dof_vec = [dof_vec;dof_inter];

end
dof_vec = unique(dof_vec);

PI = I(:,dof_vec);


