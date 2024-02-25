function out = fun_g1(x,y)

if x == 0 
    out = exp(-y);
elseif x == 1
    out = y^2+exp(-y);
elseif y == -0.25
    out = 1/16*x^2 + exp(0.25);
elseif y == 0
    out = 1;
end