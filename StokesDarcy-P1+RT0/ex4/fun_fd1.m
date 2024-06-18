function out = fun_fd1(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

% out = 2 - (x^2 - 2*x - y^2 + y + 1)/K;

out = 2 - (nu*(x^2 - 2*x - y^2 + y + 1))/K;