function out = fun_fs1(x,y,para)

nu = para.nu(x,y);

out = (cos(pi*y)*(2251799813685248*exp(x) - 3060513257434037*nu*pi^2 - 1125899906842624*nu*exp(x) + 1125899906842624*nu*pi^2*exp(x)))/1125899906842624;
