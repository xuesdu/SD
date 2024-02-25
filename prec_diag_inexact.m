function [z] = prec_diag_inexact(r, Pu, diag_invMps, diag_invMpd, omega_s, omega_d, amgParam, amgData)

    % get size
    Nu = size(Pu,1);
    Nps = length(diag_invMps);
    Npd = length(diag_invMpd);
    %Np = Nps+Npd;
    
    ru = r(1:Nu);
    rps = r(Nu+1:Nu+Nps);
    rpd = r(Nu+Nps+1:end);
    
    % preconditioning the velocity part
    %zu = Pu\ru;
    %zu = Prec_CG(Pu, ru, zeros(Nu,1), @(r)AMG_prec(r, 1, amgData, amgParam), 200, 1e-6, 1);
    zu = Prec_FGMRES(Pu, ru, zeros(Nu,1), [], @(r)AMG_prec(r, 1, amgData, amgParam), 200, 200, 1e-8, 0);
    
    % precondionting the pressure part
    zps = omega_s*(diag_invMps.*rps);
    zpd = omega_d*(diag_invMpd.*rpd);
   
    % get z
    z = [zu;zps;zpd];
       
end
