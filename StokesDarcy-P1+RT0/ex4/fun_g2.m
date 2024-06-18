function out = fun_g2(x,y)

if x == 0
    out = 2;
elseif x == 1
    out = -2/3*y^3 + 2;
elseif y == -0.25
    out = 1/96*x+2-pi*sin(pi*x);
elseif y == 0
    out = 2-pi*sin(pi*x);
end