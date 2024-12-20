function out = fun_ud2_vec(x)

out = pi*(exp(x(1,:))-x(1,:).*exp(1)).*sin(pi*x(2,:));