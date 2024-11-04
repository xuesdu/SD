function out = fun_gb2(x,y)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
global yb yt nu

if y == yb || y == yt
    out = 4*nu*exp(2*y)*exp(x) - 12*x^2*exp(y);
else
    out = 0;
end