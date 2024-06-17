function out = fun_eta2(x,y)

global K nu alpha_BJS

% out = nu*(exp(x + 2*y) + x*exp(x*y)) + (alpha_BJS*nu^(1/2)*exp(x + 2*y))/K^(1/2);


out = nu*(exp(x + 2*y) + x*exp(x*y)) + (alpha_BJS*nu*exp(x + 2*y))/K^(1/2);