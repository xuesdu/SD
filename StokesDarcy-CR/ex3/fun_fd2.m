function out = fun_fd2(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

% out = x*pi*cos(pi*y) - (x*pi*cos(pi*y))/K;

out = x*pi*cos(pi*y) - (nu*x*pi*cos(pi*y))/K;

