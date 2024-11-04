function out = fun_fd2(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);
out = (pi*sin(pi*y)*(K - 1)*(6121026514868073*x - 2251799813685248*exp(x)))/(2251799813685248*K);
% out = (cos(pi*y)*(2251799813685248*exp(x) - 6121026514868073)*(K - 1))/(2251799813685248*K);