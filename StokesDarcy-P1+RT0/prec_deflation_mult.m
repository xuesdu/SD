function [z] = prec_deflation_mult(A,PA,P,r)
% solve (P'AP)z = P'r

% preconditioner: B1 = PA
z = PA\r;           % e1 = B1*r;
r1 = r - A*z;       % r1 = r - A*e1;
% deflation --> preconditioner: B2 = P((P'AP)^{-1})P' 
z = z + P*(P'*A*P)^(-1)*P'*r1;    % e2 = e1 + B2*r1;

end

