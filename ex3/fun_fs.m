function out = fun_fs(x,y,para)
nu = para.nu(x,y);
f1 = 2*nu*pi^2*cos(pi*x)*sin(pi*y);
f2 = -pi*cos(pi*y)*(2*nu*pi*sin(pi*x) - 1);
out = [f1;f2];