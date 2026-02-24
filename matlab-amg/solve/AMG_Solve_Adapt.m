function [ x, k, err ] = AMG_Solve_Adapt(A, b, x, amgParam)
% Adaptive solve phase 
% 
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)
% 

% parameters
print_level = amgParam.print_level;
max_it = amgParam.max_it;
tol = amgParam.tol;

%setup amgData
amgParam.agg_type = 'MWM'; % this forces the aggregation on the superficial level is always HEC
[ amgData ] = AMG_Setup( A, amgParam ); %the printout of the setup phase is only for the initial AMG setup


% prepare solve 
level = 1;
err = zeros(1,1);

r = b - A*x;

err(1) = norm(r); 


if print_level > 0
    fprintf('----------------------------------------------------\n')
    fprintf('              Calling Adaptive AMG solver    \n');
    fprintf('----------------------------------------------------\n')
end

%print_level=0;
if print_level > 0
    fprintf('----------------------------------------------------\n')
    display(sprintf(' # It |  ||r||/||r0||  |     ||r||      |   Rate.  |'));
    fprintf('----------------------------------------------------\n');
    display(sprintf(' %4d |  %e  |  %e  | %f |', 0, 1.0, err(1), 0.0));
end

% main loop
solve_start = tic;

resetup = 1;
alpha = 1; % num of hierarchies using smooth error from more complicated smoothing
save_amgData = cell(1,1);
beta = 7; % num of times you repeat the complicated cycle scheme
gamma = 2; % num of recycles of all hierarchies ever generated

save_amgData{1,1} = amgData;
approx_smooth_error = zeros(length(x),1);
convRsum = 0;
w_o_restart = 0;
k = 1;
threshold = amgParam.threshold;

for k = 1:max_it
   
    % call multigrid
    x = AMG_Cycle(amgData, b, x, level, amgParam);
    
    %orthorgonalization
    x = x - x' * ones(amgData(1).N,1) / amgData(1).N;
    
    % compute residual
    r = b - A*x;
    
    % compute error
    err(k+1) = norm(r);
    
    % display
    if print_level > 0
        display(sprintf(' %4d |  %e  |  %e  | %f |', k, err(k+1)/err(1), err(k+1), err(k+1)/err(k)));
    end
    
    if ((err(end)/err(1)) < tol)
        break;
    end
    w_o_restart  = w_o_restart + 1;
    
    % re-aggregating only when the convergence becomes slower
    convRsum = err(end)/err(end-1) + convRsum;
    if convRsum / w_o_restart > threshold % for average convR
    %if err(k+1)/err(k) > threshold % this is for convR at kth step only
        if norm(b)==0
            smooth_error = x/norm(x);
            [ amgData ] = AMG_Setup_pathcover( A, smooth_error, amgParam );
        else
            if resetup < alpha + 1
                amgParam.cycle_type = 'W'; % enforce a different scheme
                approx_smooth_error = Prec_CG(A, r, approx_smooth_error, @(r)AMG_prec(r, 1, amgData, amgParam), beta, tol, 0);
                amgParam.cycle_type = 'V'; % change it back to V cycle scheme
            else
                amgParam.cycle_type = 'V'; % enforce a different scheme
                approx_smooth_error = Prec_CG(A, r, approx_smooth_error, @(r)AMG_prod_prec(r, 1, save_amgData, amgParam), gamma, tol, 0); % 2*gamma steps b/c symmetry
                amgParam.cycle_type = 'V'; % change it back to V cycle scheme
            end
            
            %orthorgonalization
            approx_smooth_error = approx_smooth_error - approx_smooth_error' * ones(amgData(1).N,1) / amgData(1).N;


            % contraction mapping constant
            CC = 1;
            x = x + CC * approx_smooth_error;
            x = x - x' * ones(amgData(1).N,1) / amgData(1).N;

            fprintf('-- resetup at iteration %3d --\n', k-1);

            [ amgData ] = AMG_Setup_pathcover( A, approx_smooth_error, amgParam );
            %[ amgData ] = AMG_Setup_MWM( A, approx_smooth_error, amgParam );

            save_amgData{resetup+1,1} = amgData; 
        end
        resetup = resetup + 1;
        w_o_restart = 0;
        convRsum = 0;
    end
    
end


solve_duration = toc(solve_start);

% print
fprintf('----------------------------------------------------\n');
if k == max_it
    fprintf('        AMG reached maximal number of iterations \n');
else
    fprintf('        AMG converged \n');
    fprintf('        Number of iterations = %d \n', k);
end
fprintf('        Relative residual = %e \n', err(end)/err(1))
fprintf('----------------------------------------------------\n');
fprintf('        AMG solve costs %f seconds\n', solve_duration);
fprintf('----------------------------------------------------\n');





end
