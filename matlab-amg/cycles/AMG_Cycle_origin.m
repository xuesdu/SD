function [x] = AMG_Cycle(amgData, b, x, level, amgParam)
% Multigrid cycle
%
% @ Xiaozhe Hu, Tufts University
    
% parameters
max_level = amgData(1).max_level;

cycle_type = amgParam.cycle_type;

n_presmooth = amgParam.n_presmooth;
n_postsmooth = amgParam.n_postsmooth;

smooth_type = amgParam.smooth_type;
omega = amgParam.smooth_omega;
randomized = 1;

ILU_level = amgParam.ILU_level;

% coarsest level
if (level == max_level) 
    % x = NA( amgData(level).A, b, amgParam.coarse_it);   
    x = amgData(level).A \ b;
    %n_c = size(b,1);
    %x = ( amgData(level).A + 1.0e-6*speye(n_c,n_c) ) \ b;
    %x = ( amgData(level).A + 1.0e-10*speye(length(b),length(b)) ) \ b;
    %x = jacobi(amgData(level).A, b, x, amgData(level).Dinv, 1);
    %x = forward_gs(amgData(level).A, b, x, amgData(level).DL, 1);
    %x = backward_gs(amgData(level).A, b, x, amgData(level).DL, 1);

else
    
    % presmoothing
    if level <= ILU_level
        x = ILU_smoother(amgData(level).A, b, x, amgData(level).IL, amgData(level).IU, n_presmooth);
    else
        switch smooth_type
            case 'Rd',
                omega = 1/amgData(level).infNorm;
                x = richardson(amgData(level).A, b, x, omega, n_presmooth);
            case 'Jcb',
                x = jacobi(amgData(level).A, b, x, amgData(level).Dinv, n_presmooth);
            case 'K1',
                x = kaczmarz_AAt(amgData(level).A, b, x, n_presmooth, randomized);
            case 'K2',
                x = kaczmarz_AtA(amgData(level).A, b, x, n_presmooth, randomized);
            otherwise,
                x = forward_gs(amgData(level).A, b, x, amgData(level).DL, n_presmooth);
        end
    end
    
    % compute residual
    r = b - amgData(level).A*x;
    
    % restriction
    r_c = amgData(level).R*r;
                
    % coarse grid correction
    [n_c, m_c] = size(r_c);
    e_c = zeros(n_c,m_c);
    
    switch cycle_type
    
        case 'V',  
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            
        case 'W',
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            
        case 'nV',  
            for i = 1:amgParam.coarse_it
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            end
            
        case 'N',
            if (level == max_level-1)
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            else
                L = 1 ;  
                %mu = 0;
                mu = ( 1 -  0.227928  );%/( 1 + 0.852942 ) ;    %
                %      32         64         128        256        512     1024
                %  0.731077   0.809767    0.833944   0.852942   0.834982   0.845074  MIS %Vcycle  amgParam.max_level = 3;amgParam.coarse_it = 1;
                %  0.515110   0.625238    0.627858   0.634707               HEC
                %  0.214256   0.224779    0.226951   0.227928         MWM
                alpha = 1 / L ;
                beta = (sqrt(L) - sqrt(mu))/(sqrt(L) + sqrt(mu)) ; 
                e_0 = zeros(n_c,m_c);
                
                p = AMG_Cycle(amgData, r_c, e_0, level+1, amgParam); 
                e_c_old = e_c;
                %e_c = alpha*p;
                %e_c = zeros(n_c,m_c);
                e_c = (r_c'*p)*p / (p'*amgData(level+1).A*p); 
                %e_c = ( 1 + beta ) * (e_c_old+alpha*p);

                for i = 2: amgParam.coarse_it
                    
                    q = AMG_Cycle(amgData, r_c- amgData(level+1).A*e_c, e_0, level+1, amgParam);
                    
                    e_c_new = ( 1 + beta ) * (e_c+alpha*q) - beta * (e_c_old+alpha*p); 
                    
                    e_c_old = e_c;
                    e_c = e_c_new;
                    p = q;
                    
                end
            end
            
        case 'H'
             L = 1 ;  
             %mu = 0.1;
             mu = ( 1 -   0.227928);%/( 1 + 0.731077  ) ; 
             alpha = 4 / ((sqrt(L)+sqrt(mu))^2); 
             beta = ( max( abs( 1 - sqrt(alpha * mu) ), abs( 1 - sqrt(alpha*L) ) ) )^2;
             e_0 = zeros(n_c,m_c);
             p = AMG_Cycle(amgData, r_c, e_0, level+1, amgParam);
             e_c_old = e_c;
             %e_c = alpha*p;
              e_c = (r_c'*p)*p / (p'*amgData(level+1).A*p); 
             for i = 2: amgParam.coarse_it
                 
                    q = AMG_Cycle(amgData, r_c- amgData(level+1).A*e_c, e_0, level+1, amgParam);
                    
                    e_c_new =  e_c + alpha * q  + beta * ( e_c - e_c_old); 
                    
                    e_c_old = e_c;
                    e_c = e_c_new; 
             end
           
        case 'Cold'   %% Chebychev polynomial used in semi iterative method
             C = zeros(amgParam.coarse_it,1); 
             C(1) = 1 ; 
             rho =  0.515110 ; 
             C(2) = 1/rho ;
             for i = 3:amgParam.coarse_it
                 C(i) = 2 * C(i-1) / rho - C(i-2);
             end
             
             e_0 = zeros(n_c,m_c);
             p = AMG_Cycle(amgData, r_c, e_0, level+1, amgParam);
             e_c_old = e_c;
             e_c = p;
             %e_c = 2 * ( C(1) / (rho * C(2)) ) * p ;
             for i = 2: amgParam.coarse_it
                 
                    q = AMG_Cycle(amgData, r_c - amgData(level+1).A*e_c, e_0, level+1, amgParam);
                    
                    e_c_new = 2 * ( C(i-1) / (rho * C(i)) ) * ( e_c +  q - e_c_old) +  e_c_old ; 
                    
                    e_c_old = e_c;
                    e_c = e_c_new; 
             end          
        
        case 'C' %AMLI-cycle using scaled and shifted Chebyshev polynomial
            if (level == max_level-1)
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            else
                L = 1; mu = 0.1; e_0 = e_c;
                switch amgParam.coarse_it
                    case 2
                        w=2/(L+mu);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);

                    case 3
                        w = 4/(L+3*mu);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        w = 1/L;
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);

                    case 4
                        w = 2*sqrt(2)/((sqrt(2)-1)*L + (sqrt(2)+1)*mu);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        w = 2*sqrt(2)/((sqrt(2)+1)*L + (sqrt(2)-1)*mu);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);

                    case 5
                        w = 2/((L+mu)-cos(pi/5)*(L-mu));
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        w = 2/((L+mu)-cos(3*pi/5)*(L-mu));
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        w = 1/L;
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        
                    otherwise
                        for i=1:2:(amgParam.coarse_it-1)
                            w = 2/((L+mu)-cos(i*pi/amgParam.coarse_it)*(L-mu));
                            e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        end
   
                        if ((i+2)==amgParam.coarse_it)
                            w = 1/L;
                            e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        end
                end

            end

        case 'M' %AMLI-cycle using momemtum acceleration polynomial 
            if (level == max_level-1)
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            else
                a = 1.99; e_0 = e_c;
                switch amgParam.coarse_it
                    case 2
                        L = (2+a)^2/(8*a);
                        w = 2/L;
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);
                        w = a/L;
                        e_c = e_c + w*AMG_Cycle(amgData, r_c-amgData(level+1).A*e_c, e_0, level+1, amgParam);

                end
            end
            
        case 'NAV' 
            if (level == max_level-1)
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            else
                e_0 = zeros(n_c,m_c);
                e_c =  NA_prec(amgData(level+1).A, r_c,  amgParam.coarse_it,...
                   @(r)AMG_Cycle(amgData, r, e_0, level+1, amgParam));
%                 e_c =  NA_CG(amgData(level+1).A, r_c,  amgParam.coarse_it,...
%                    @(r)AMG_Cycle(amgData, r, e_0, level+1, amgParam));
            end

        case 'K',
            if (level == max_level-1)
                e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
            else 
                e_0 = zeros(n_c,m_c);
                switch amgParam.coarse_krylov_type
                    case 'GCG',                    
                        e_c = gcg(amgData(level+1).A, r_c, 1e-12, amgParam.coarse_it, ...
                            @(r)AMG_Cycle(amgData, r, e_0, level+1, amgParam), 0);
                    otherwise,
                        e_c = Prec_FGMRES(amgData(level+1).A, r_c, e_c, [], ...
                            @(r)AMG_Cycle(amgData, r, e_0, level+1, amgParam),...
                            amgParam.coarse_it, amgParam.coarse_it, 1e-12, 0);
                end
            end
            
        otherwise,
            display('Wrong cycle type!! run V-cycle!!')
            e_c = AMG_Cycle(amgData, r_c, e_c, level+1, amgParam);
    end
            
    % scaling
    alpha = 1.0;
    %alpha = (r_c'*e_c)/(e_c'*amgData(level+1).A*e_c);
    
    % prolongation
    x = x + alpha*amgData(level).P*e_c;
        
    % postsmoothing
    if level <= ILU_level
        x = ILU_smoother(amgData(level).A, b, x, amgData(level).IL, amgData(level).IU, n_postsmooth);
    else
        switch smooth_type
            case 'Rd',
                omega = 1/amgData(level).infNorm;
                x = richardson(amgData(level).A, b, x, omega, n_postsmooth);
            case 'Jcb',
                x = jacobi(amgData(level).A, b, x, amgData(level).Dinv, n_postsmooth);
            case 'K1',
                x = kaczmarz_AAt(amgData(level).A, b, x, n_postsmooth, randomized);
            case 'K2',
                x = kaczmarz_AtA(amgData(level).A, b, x, n_postsmooth, randomized);
            otherwise,
                x = backward_gs(amgData(level).A, b, x, amgData(level).DU, n_postsmooth);
        end
    end
        
end

end
