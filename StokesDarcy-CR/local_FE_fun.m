function out = local_FE_fun(x,y,uh_local,vertices,basis_type,basis_der_x,basis_der_y)
% uh = sum_{j=1}^{N+1} u_j*phi_j
out = 0;

number_of_local_basis_trial = length(uh_local);

for k = 1:number_of_local_basis_trial
    out = out + uh_local(k)*FE_basis_fun_local(x,y,vertices,basis_type,k,basis_der_x,basis_der_y);
end
