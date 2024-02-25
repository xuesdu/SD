function out = fun_fs2(x,y,para)

nu = para.nu(x,y);
out = -(exp(x)*sin(pi*y)*(nu*pi^2 - nu + 2*pi^2))/pi;