function out = fun_fd1(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

% out = exp(x*y)/K + 16*y^3;

out = 16*y^3 + (nu*exp(x*y))/K;
