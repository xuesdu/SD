function out = fun_fd2(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

%out = ((2*y - 1)*(x - 1))/K + 2;

out = (nu*(2*y - 1)*(x - 1))/K + 2;