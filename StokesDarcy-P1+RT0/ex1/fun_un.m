function out = fun_un(x,y)

global xr yb yt xbar

if x == xr
    out = -(exp(x)-exp(1)).*cos(pi*y);
elseif x == xbar
    out = (exp(x)-exp(1)).*cos(pi*y);
elseif y == yb
    out = -pi*(exp(x)-x*exp(1)).*sin(pi*y);
elseif y == yt
    out = pi*(exp(x)-x*exp(1)).*sin(pi*y);
else
    fprintf('error');
end