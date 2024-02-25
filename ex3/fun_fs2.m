function out = fun_fs2(x,y,para)

nu = para.nu(x,y);
out = -pi*cos(pi*y)*(2*nu*pi*sin(pi*x) - 1);