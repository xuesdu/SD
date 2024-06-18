function out = fun_fs(x,y,para)
nu = para.nu(x,y);
f1 = 24*x*exp(y) - 2*nu*(exp(x + 2*y) + (x^2*exp(x*y))/2 + y^2*exp(x*y));
f2 = 12*x^2*exp(y) - 2*nu*(exp(x*y)/2 + (9*exp(2*y)*exp(x))/2 + (x*y*exp(x*y))/2);
out = [f1;f2];