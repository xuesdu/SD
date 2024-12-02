function out = fun_fd1(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

out = cos(pi*y)*(exp(x) - 6121026514868073/2251799813685248) - (cos(pi*y)*(exp(x) - 6121026514868073/2251799813685248) - 1)/K;