function out = fun_fs(x,y,para)
nu = para.nu(x,y);
f1 = 2*exp(x)*cos(pi*y) - nu*(exp(x)*cos(pi*y) - pi^2*cos(pi*y)*(exp(x) - 3060513257434037/1125899906842624));
f2 = - nu*(pi*exp(x)*sin(pi*y) - (exp(x)*sin(pi*y))/pi) - 2*pi*exp(x)*sin(pi*y);
out = [f1;f2];