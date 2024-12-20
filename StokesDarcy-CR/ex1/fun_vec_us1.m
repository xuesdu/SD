function out = fun_vec_us1(x)

out = (exp(x(1,:))-exp(1)).*cos(pi*x(2,:));