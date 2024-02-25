function out = fun_fd(x,y,para)

K = para.K(x,y);

f1 = sin(pi*y) - sin(pi*y)/K;
f2 = x*pi*cos(pi*y) - (x*pi*cos(pi*y))/K;
out = [f1;f2];