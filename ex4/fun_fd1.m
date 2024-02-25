function out = fun_fd1(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

out = 3*x^2*y^3;
