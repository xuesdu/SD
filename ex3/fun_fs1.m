function out = fun_fs1(x,y,para)

nu = para.nu(x,y);

out = 2*nu*pi^2*cos(pi*x)*sin(pi*y);