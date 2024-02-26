function [z] = prec_diag_exact(r, Pu, diag_invMps, diag_invMpd, omega_s, omega_d)

    % get size
    Nu = size(Pu,1);
    Nps = length(diag_invMps);
    Npd = length(diag_invMpd);
    %Np = Nps+Npd;
    
    ru = r(1:Nu);
    rps = r(Nu+1:Nu+Nps);
    rpd = r(Nu+Nps+1:end);
    
    % preconditioning 
    %zL = Prec_FGMRES(ML, rL, zeros(NL,1), [], @(r)AMG_prec(r, 1, amgData_ML, amgParam), 100, 100, 1e-2, 0);
    zu = Pu\ru;
    %zP = Prec_FGMRES( S, rP, zeros(NP,1), [], @(r)AMG_prec(r, 1, amgData_S,  amgParam), 100, 100, 1e-2, 0);
    zps = omega_s*(diag_invMps.*rps);
    zpd = omega_d*(diag_invMpd.*rpd);
   
    % get z
    z = [zu;zps;zpd];
       
end
