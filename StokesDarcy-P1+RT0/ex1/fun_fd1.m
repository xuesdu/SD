function out = fun_fd1(x,y,para)

K = para.K(x,y);
nu = para.nu(x,y);

out = (cos(pi*y)*(K - nu)*(1125899906842624*exp(x) - 3060513257434037))/(1125899906842624*K);