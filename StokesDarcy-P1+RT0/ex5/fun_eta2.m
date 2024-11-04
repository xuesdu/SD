function out = fun_eta2(x,y)

global nu K alpha_BJS

out = - (nu*sin(pi*y)*(1125899906842624*exp(x) - 3060513257434037*pi^2 + 1125899906842624*pi^2*exp(x)))/(1125899906842624*pi) - (alpha_BJS*nu^(1/2)*exp(x)*sin(pi*y))/(K^(1/2)*pi);

%out = - nu*((exp(x)*sin(pi*y))/pi + pi*sin(pi*y)*(exp(x) - 3060513257434037/1125899906842624)) - (alpha_BJS*nu*exp(x)*sin(pi*y))/(K^(1/2)*pi);