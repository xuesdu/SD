function out = fun_eta2(x,y)

global nu K alpha_BJS

out = - nu*((exp(x)*sin(pi*y))/pi + pi*sin(pi*y)*(exp(x) - 3060513257434037/1125899906842624)) - (alpha_BJS*nu*exp(x)*sin(pi*y))/(K^(1/2)*pi);
