function out = fun_gb2(x,y,para)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
nu = para.nu(x,y);

out = 2 - 2*y - 2*x;