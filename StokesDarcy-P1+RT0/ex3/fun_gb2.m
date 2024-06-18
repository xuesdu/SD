function out = fun_gb2(x,y,para)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
nu = para.nu(x,y);

out = 2*nu*pi*sin(pi*x)*sin(pi*y) - sin(pi*y);