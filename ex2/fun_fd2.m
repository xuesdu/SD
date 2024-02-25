function out = fun_fd2(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

% out = (y^2 + 3*x)/K + 48*x*y^2;

out = 48*x*y^2 + (nu*(y^2 + 3*x))/K;
