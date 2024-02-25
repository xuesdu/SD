function out = fun_eta1(x,y)

global nu
out = 12*x^2*exp(y) - 16*x*y^3 - 2*nu*y*exp(x*y) + 5312313071119285/1125899906842624;
