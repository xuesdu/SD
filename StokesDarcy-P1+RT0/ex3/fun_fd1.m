function out = fun_fd1(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

% out = sin(pi*y) - sin(pi*y)/K;


 out = sin(pi*y) - (nu*sin(pi*y))/K;