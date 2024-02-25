function out = fun_fd2(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

out = 3*x^3*y^2;
