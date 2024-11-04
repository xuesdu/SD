function out = fun_gb2(x,y)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
global yb yt nu

if y == yb || y == yt
    out = 2*nu*pi*sin(pi*x)*sin(pi*y) - sin(pi*y);
else
    out = 0;
end