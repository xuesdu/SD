function out = fun_gb2(x,y)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
global yb yt nu

if y == yb || y == yt
    out = - 2*exp(x)*cos(pi*y) - 2*nu*exp(x)*cos(pi*y);
else
    out = 0;
end