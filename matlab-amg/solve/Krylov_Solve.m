function [ x, iter, residual, alpha_all, beta_all] = Krylov_Solve( A, b, x, iterParam, amgParam, amgData)
% Solve phase using Krylov method
%
% @ Xiaozhe Hu, Tufts University

max_it = iterParam.max_it;
tol = iterParam.tol;
print_level = iterParam.print_level;
restart = iterParam.restart;

alpha_all = [];
beta_all = [];

if norm(b,'fro') < 1e-12
    x = zeros(length(b),1);
    if print_level > 0
        fprintf('----------------------------------------------------\n')
        fprintf('     zero right hand side, set solution to zero    \n');
        fprintf('----------------------------------------------------\n')
    end
    iter = 0;
    residual = 0;
    return;
end

if print_level > 0
    fprintf('----------------------------------------------------\n');
end

solve_start = tic;

switch iterParam.solver_type
    
    case 'CG',
        
        if print_level > 0
            fprintf('        Calling Preconditioned CG solver    \n');
        end
        
        switch iterParam.prec_type
            
            case 'Jacobi'
                n = size(A,1);
                D = spdiags(diag(A), 0, n, n);
                [x, iter, residual, alpha_all, beta_all] = Prec_CG(A, b, x, D, max_it, tol, print_level);
            
            case 'AMG'
                
                [x, iter, residual, alpha_all, beta_all] = Prec_CG(A, b, x, @(r)AMG_prec(r, 1, amgData, amgParam), max_it, tol, print_level);
                %[x, ~, ~, iter, residual]  = pcg(A, b,tol, max_it, @(r)AMG_prec(r, 1, amgData, amgParam));
                                
            otherwise
                
                if print_level > 0
                    display('  No preconditioner sepecified, run CG')
                end
                [x, iter, residual, alpha_all, beta_all] = Prec_CG(A, b, x, [], max_it, tol, print_level);
                %[x, ~, ~, iter, residual]  = pcg(A, b,tol, max_it);
                
        end
        
    case 'FGMRES',
        
        if print_level > 0
             fprintf('       Calling Preconditioned GMRES solver    \n');
        end
        
        switch iterParam.prec_type
            
            case 'AMG'
                
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], @(r)AMG_prec(r, 1, amgData, amgParam), max_it, restart, tol, print_level);
                
            otherwise
                
                if print_level > 0
                    display('       No preconditioner sepecified, run FGMRES');
                end
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], [], max_it, restart, tol, print_level);
        end   
        
    otherwise

        display('Wrong solver type!!')
end

solve_duration = toc(solve_start);
  
if print_level > 0
    fprintf('----------------------------------------------------\n');

    fprintf('----------------------------------------------------\n');
end

if print_level > -1 
    if iter == max_it
            fprintf('        Krylov method reached maximal number of iterations \n');
    else
            fprintf('        Krylov method converged \n');
            fprintf('        Number of iterations = %d \n', iter);
            fprintf('----------------------------------------------------\n');
            fprintf('        Krylov method costs %f seconds\n', solve_duration);
            fprintf('----------------------------------------------------\n');
    end
end

if print_level > 0
    fprintf('        Relative residual = %e \n', residual(iter)/residual(1))
end

end

