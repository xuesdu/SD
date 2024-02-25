function out = fun_fd(x,y,para)

K = para.K(x,y);

f1 = exp(x*y)/K + 16*y^3;
f2 = (y^2 + 3*x)/K + 48*x*y^2;
out = [f1;f2];