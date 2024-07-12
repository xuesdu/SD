function out = fun_eta2(x,y)

global nu K alpha_BJS

out = nu*(2*x + 2*y - 3) + (alpha_BJS*nu^(1/2)*(x - 1)^2)/K^(1/2);


% out = nu*(2*x + 2*y - 3) + (alpha_BJS*nu*(x - 1)^2)/K^(1/2);