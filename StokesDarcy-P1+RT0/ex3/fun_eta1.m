function out = fun_eta1(x,y)
global nu
out = sin(pi*y)*(2*nu*pi*sin(pi*x) - x + 1);