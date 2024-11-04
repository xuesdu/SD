function out = fun_gb1(x,y)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
global yb yt nu

if y == yb || y == yt
    out = nu*(exp(2*y)*exp(x) + x*exp(x*y));
else
    out = 0;
end