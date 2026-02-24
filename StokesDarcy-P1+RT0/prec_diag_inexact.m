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
    % [zu, iter_inner, residual] = Prec_FGMRES(Pu, ru, zeros(Nu,1), [], @(r)AMG_prec(r, 1, amgData, amgParam), amgParam.max_it, amgParam.max_it, amgParam.tol, 0);
    
    [zu, iter_inner, residual] = AMG_Solve(amgData, ru, zeros(Nu,1), amgParam);
    fprintf('%d\n', iter_inner);
            
    
    % if (amgParam.print_level == 0)
    % 
    %     if (iter_inner == amgParam.max_it)
    %         fprintf('----------------------------------------------------\n');
    %         fprintf('   UA-AMG reached maximal number of iterations \n');
    %         fprintf('----------------------------------------------------\n');
    %     else
    %         fprintf('----------------------------------------------------\n');
    %         fprintf('UA-AMG converged at iteration %d with relative residual %e\n', iter_inner, residual(iter_inner+1)/residual(1));
    %         fprintf('----------------------------------------------------\n');
    % 
    %     end
    % end
    % precondionting the pressure part
    zps = omega_s*(diag_invMps.*rps);
    zpd = omega_d*(diag_invMpd.*rpd);
   
    % get z
    z = [zu;zps;zpd];
       
end
