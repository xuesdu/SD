function out = fun_fs2(x,y,para)

nu = para.nu(x,y);

out = 12*x^2*exp(y) - 2*nu*(exp(x*y)/2 + (9*exp(2*y)*exp(x))/2 + (x*y*exp(x*y))/2);
