function out = fun_eta1(x,y)
global nu

out = (cos(pi*y)*(3060513257434037*x + 1125899906842624*exp(x) - 2251799813685248*nu*exp(x)))/1125899906842624;
