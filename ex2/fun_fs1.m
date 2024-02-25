function out = fun_fs1(x,y,para)

nu = para.nu(x,y);

out = 24*x*exp(y) - 2*nu*(exp(x + 2*y) + (x^2*exp(x*y))/2 + y^2*exp(x*y));