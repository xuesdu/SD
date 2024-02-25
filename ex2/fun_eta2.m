function out = fun_eta2(x,y)

global K

out = - exp(x + 2*y) - K^(1/2)*(exp(x + 2*y) + x*exp(x*y));