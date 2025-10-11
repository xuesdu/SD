function [z] = prec_deflation_add_copy(lambda,A,B1,P_kernel,PA_kernel,r)
% solve (P'AP)z = P'r

% preconditioner: B1 = PA
z1 = B1\r;   % z = B1*r;
% deflation --> preconditioner: B2 = P((P'AP)^{-1})P' 
z2 = PA_kernel*lambda^(-1)*(P_kernel'*A*P_kernel)^(-1)*PA_kernel'*r; % z = B2*r
z = z1 + z2;

end

