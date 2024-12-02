function out = fun_eta1(x,y)
global nu

out = (6121026514868073*x*cos(pi*y))/2251799813685248 - 2*nu + exp(x)*cos(pi*y);
