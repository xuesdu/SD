function out = fun_gbp(x,y,para)
% p = gbp; Neumann boundary condition
nu = para.nu(x,y);

out = (exp(x)-x*exp(1))*cos(pi*y);