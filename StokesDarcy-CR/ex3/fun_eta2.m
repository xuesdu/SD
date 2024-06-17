function out = fun_eta2(x,y)

global alpha_BJS nu K

% out = -(alpha_BJS*nu^(1/2)*cos(pi*y)*sin(pi*x))/K^(1/2);


out = -(alpha_BJS*nu*cos(pi*y)*sin(pi*x))/K^(1/2);

