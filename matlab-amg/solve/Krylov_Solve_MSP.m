function [ x, iter, residual ] = Krylov_Solve( A, H, M, b, x, iterParam, amgParam, amgData, P_mat)
% Solve phase using Krylov method
%
% @ Xiaozhe Hu, Tufts University

max_it = iterParam.max_it;
tol = iterParam.tol;
print_level = iterParam.print_level;
restart = iterParam.restart;

if norm(b) < 1e-12
    x = zeros(length(b),1);
    if print_level > 0
        fprintf('----------------------------------------------------\n')
        fprintf('     zero right hand side, set solution to zero    \n');
        fprintf('----------------------------------------------------\n')
    end
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
                [x, iter, residual] = Prec_CG(A, b, x, D, max_it, tol, print_level);
                
            case 'SGS'
                n = size(A,1);
                D_inv = spdiags(diag(A).^(-1), 0, n, n);
                DL = tril(A);
                DU = triu(A);
                %D_SGS = DL * D_inv * DU;
                %[x, iter, residual] = Prec_CG(A, b, x, D_SGS, max_it, tol, print_level);
                [x, iter, residual] = Prec_CG(A, b, x, @(r)SGS_prec(r, DL, DU, D_inv), max_it, tol, print_level);
                
            case 'ILU'
                %n = size(A,1);
                ILU_setup.type = amgParam.ILU_type;  % nofill, crout, ilutp
                ILU_setup.droptol = amgParam.droptol;
                ILU_setup.milu = amgParam.ILU_milu;  % row, col, off
                ILU_setup.udiag = amgParam.ILU_udiag;
                ILU_setup.thresh = amgParam.ILU_thresh;
                
                [IL, IU] = ilu(A, ILU_setup);
                %[x, iter, residual] = Prec_CG(A, b, x, IL*IU, max_it, tol, print_level);
                [x, iter, residual] = Prec_CG(A, b, x, @(r)ILU_prec(r, IL, IU), max_it, tol, print_level);
                
            case 'AMG'
                %n = size(A,1);
                %[x, iter, residual] = Prec_CG(A, b, x, H+1e-6*eye(n), max_it, tol, print_level);
                [x, iter, residual] = Prec_CG(A, b, x, @(r)AMG_prec(r, 1, amgData, amgParam), max_it, tol, print_level);
                %[x, ~, ~, iter, residual]  = pcg(A, b,tol, max_it, @(r)AMG_prec(r, 1, amgData, amgParam));
                   
            case 'MSP'
                n = size(M,1);
                [x, iter, residual] = Prec_CG(A, b, x, @(r)MSP_prec(r, M+1e-6*speye(n), H+1e-6*speye(n), P_mat, 5), max_it, tol, print_level);
%                 [x, iter, residual] = Prec_CG(A, b, x, H+1e-6*eye(n), max_it, tol, print_level);
            case 'blkGE'
                n = length(x);
                [x, iter, residual] = Prec_CG(A, b, x, @(r)blkGE_prec(r, H,  n, iterParam.layer_ind, iterParam.max_level, 1), max_it, tol, print_level);
%                 n = size(A,1);
%                 [x, iter, residual] = Prec_CG(A, b, x, H+1e-6*eye(n), max_it, tol, print_level);
            otherwise
                
                if print_level > 0
                    display('  No preconditioner sepecified, run CG')
                end
                %[x, iter, residual] = Prec_CG(A, b, x, [], max_it, tol, print_level);
                try
                    n = size(M,1);
                    %[x, iter, residual] = Prec_CG(A, b, x, M+1e-6*speye(n), max_it, tol, print_level);
                    [x, iter, residual] = Prec_CG(A, b, x, @(r)MSP_prec(r, M+1e-6*speye(n), [], P_mat, []), max_it, tol, print_level);
                    %[x, iter, residual] = Prec_CG(A, b, x, @(r)MSP_prec(r, M, [], P_mat, []), max_it, tol, print_level);


                catch
                    n = size(H,1);
                    %[x, iter, residual] = Prec_CG(A, b, x, @(r)MSP_prec(r, [],H, P_mat, []), max_it, tol, print_level);
                    [x, iter, residual] = Prec_CG(A, b, x, @(r)MSP_prec(r, [],H+1e-6*speye(n), P_mat, []), max_it, tol, print_level);
                    %[x, iter, residual] = Prec_CG(A, b, x, H, max_it, tol, print_level);
                end
                    %[x, ~, ~, iter, residual]  = pcg(A, b,tol, max_it);
                
        end
        
    case 'FGMRES',
        
        if print_level > 0
             fprintf('       Calling Preconditioned GMRES solver    \n');
        end
        
        switch iterParam.prec_type
            
            case 'AMG'
                
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], @(r)AMG_prec(r, 1, amgData, amgParam), max_it, restart, tol, print_level);
                
            case 'MSP'
                n = size(M,1);
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], @(r)MSP_prec(r, M+1e-6*speye(n), H+1e-6*speye(n), P_mat, 5), max_it, restart, tol, print_level);

            otherwise
                
                if print_level > 0
                    display('       No preconditioner sepecified, run FGMRES');
                end
                %n = size(A,1);
                %[x, iter, residual] = Prec_FGMRES(A, b, x, [], H, max_it, restart, tol, print_level);
                [x, iter, residual] = Prec_FGMRES(A, b, x, [], [], max_it, restart, tol, print_level);
        end   
    case 'GCG',
        
        if print_level > 0
             fprintf('       Calling Preconditioned GCG solver    \n');
        end
        switch iterParam.prec_type
            
            case 'AMG'
                
                [x, iter, residual] = gcg2(A, b, x, tol, max_it, @(r)AMG_prec(r, 1, amgData, amgParam), print_level);
                
            case 'MSP'
                n = size(M,1);
                [x, iter, residual] = gcg2(A, b, x, tol, max_it, @(r)MSP_prec(r, M+1e-6*speye(n), H+1e-6*speye(n), P_mat, 10),  print_level);

                %[x, iter, residual] = Prec_FGMRES(A, b, x, [], @(r)MSP_prec(r, M, H, 1e-3), max_it, restart, tol, print_level);

            otherwise
                
                if print_level > 0
                    display('       No preconditioner sepecified, run GCG');
                end
                %n = size(A,1);
                %[x, iter, residual] = gcg2(A, b, x, tol, max_it, H, print_level);
                [x, iter, residual] = gcg2(A, b, x, tol, max_it, [], print_level);
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

