function out = fun_gb2(x,y,para)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
nu = para.nu(x,y);

out = 4*nu*exp(2*y)*exp(x) - 12*x^2*exp(y);