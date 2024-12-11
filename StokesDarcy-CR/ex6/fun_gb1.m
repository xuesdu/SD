function out = fun_gb1(x,y,para)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
nu = para.nu(x,y);

out = -nu*((exp(x)*sin(pi*y))/pi - 2);