function out = fun_gb1(x,y)
% sigma_s\cdot n_s = [gb1; gb2]; Neumann boundary condition
global yb yt nu

if y == yb || y == yt
    out = -nu*((exp(x)*sin(pi*y))/pi + pi*sin(pi*y)*(exp(x) - 3060513257434037/1125899906842624));
else
    out = 0;
end