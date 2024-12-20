function out = fun_ud1_vec(x)

out = -(exp(x(1,:))-exp(1)).*cos(pi*x(2,:));