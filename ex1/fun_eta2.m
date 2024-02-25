function out = fun_eta2(x,y)

global K alpha

% out = K^(1/2)*((exp(x)*sin(pi*y))/pi + pi*sin(pi*y)*(exp(x) - 3060513257434037/1125899906842624)) + (exp(x)*sin(pi*y))/pi;


out = (K^(1/2)*((exp(x)*sin(pi*y))/pi + pi*sin(pi*y)*(exp(x) - 3060513257434037/1125899906842624)))/alpha + (exp(x)*sin(pi*y))/pi;