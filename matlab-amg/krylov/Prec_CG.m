function [ u, k, residual, alpha_all, beta_all] = Prec_CG(A, f, u, M, maxit, tol, print_level)
% Preconditioned Conjugate Gradient Method 
%
% @ Xiaozhe Hu, Tufts University

%-------------------
% Preparation 
%-------------------
% size of the problem
N = size(f,1);

r = zeros(N,1);
%p = zeros(N,1);
z = zeros(N,1);
Ap = zeros(N,1);
residual = zeros(maxit+1,1);

% store all alpha and beta for estimating the condition numbers
alpha_all = zeros(maxit,1);
beta_all = zeros(maxit,1);
%rz = 0.0;

if isa(A, 'double')
    r = f - A*u;
elseif isa(A, 'function_handle')
    r = f - A(u);
else
    error('A is neither a matrix or a function handle!!!');
end % end if
normr = norm(r);
residual(1) = normr;

if (normr < 10^-10)
    return;
end

if (print_level > 1)
    fprintf('----------------------------------------------------\n');
    fprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |\n');
    fprintf('----------------------------------------------------\n');
    fprintf(' %4d |  %e  |  %e  | %f |\n', 0, 1.0, residual(1), 0.0);
end

% Preconditioning: z = M\r
if isempty(M)
    z = r;
elseif isa(M, 'double')
    z = M\r;
elseif isa(M, 'function_handle')
    z = M(r);
else
    error('Preconditoner M is invalid!!!');
end % end if

p = z;

rz = r'*z;

%-------------------
% Main loop 
%-------------------
for k = 1:maxit
    
    % Ap
    if isa(A, 'double')
        Ap = A*p;
    elseif isa(A, 'function_handle')
        Ap = A(p);
    else
        error('A is neither a matrix or a function handle!!!');
    end % end if
    
    % alpha = (r,z)/(Ap,p)
    alpha = rz/(Ap'*p);
    
    % store alpha
    alpha_all(k) = alpha;
    
    % u = u + alpha*p
    u = u + alpha*p;
    
    % r = r - alpha*Ap
    r = r - alpha*Ap;
    residual(k+1) = norm(r);
 
    % display
    if (print_level > 1)
        fprintf(' %4d |  %e  |  %e  | %f |\n', k, residual(k+1)/residual(1), residual(k+1), residual(k+1)/residual(k));
    end
        
    if ((residual(k+1)/residual(1)) < tol)
        break;
    end
    
    % Preconditioning: z = M\r
    if isempty(M)
        z = r;
    elseif isa(M, 'double')
        z = M\r;
    elseif isa(M, 'function_handle')
        z = M(r);
    else
        error('Preconditoner M is invalid!!!');
    end % end if
    
    % (r_new, z_new)
    rz_new = r'*z;
    
    % beta = (r_new, z_new)/(r,z)
    beta = rz_new/rz;
    rz = rz_new;
    
    % store beta
    beta_all(k) = beta;
    
    % p = z + beta*p
    p = z + beta*p;
    
end %for

% cut residual
iter = k;
residual = residual(1:iter+1);

% cut stored alpha and beta
alpha_all = alpha_all(1:iter);
beta_all = beta_all(1:iter-1);

% print 
if (print_level > 0)
   
    if (iter == maxit)
            fprintf('----------------------------------------------------\n');
            fprintf('   CG reached maximal number of iterations \n');
            fprintf('----------------------------------------------------\n');
    else
            fprintf('----------------------------------------------------\n');
            fprintf('CG converged at iteration %d with relative residual %e\n', iter, residual(iter+1));
            fprintf('----------------------------------------------------\n');

    end
    
end

end

