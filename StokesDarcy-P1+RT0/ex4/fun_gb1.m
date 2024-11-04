function out = fun_gb1(x,y)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
global yb yt nu

if y == yb || y == yt
    out = nu*(2*x + 2*y - 3);
else
    out = 0;
end