function [z] = prec_deflation_add(A, PA, P, r)
% solve (P'AP)z = P'r

% preconditioner: B1 = PA
z1 = PA\r;            % z1 = z0 + B1*r;
% deflation --> preconditioner: B2 = P((P'AP)^{-1})P' 
z2 = P*(P'*A*P)^(-1)*P'*r; % z = B2*r
z = z1 + z2;

end

