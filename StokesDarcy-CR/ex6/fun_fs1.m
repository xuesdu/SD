function out = fun_fs1(x,y,para)

nu = para.nu(x,y);

out = exp(x)*cos(pi*y)*(nu + 2);
