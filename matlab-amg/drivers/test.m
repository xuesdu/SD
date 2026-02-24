
% test AMG
% 
% @ Xiaozhe Hu, Tufts University

%---------------------
% generate matrices
%---------------------
tic;
% grid-like graphs
% N = 64^2;
% L = assembleGraphLaplace(sqrt(N));
%N = 5;
%L = gallery('tridiag', N, -1,2,-1); 
%h = 1/(N+1);
%L = (1/h)*L;
%L = AN(1:1296,1:1296);
L = A;

% E-R random graphs
%N = 4000;
%p = 2*log(N)/N;
%[Adj, G] = Erdos_Reyni_Random_Graph( N, p );
%L = laplacian(G); % N^2

% B = incidence(G);
% Bt = B';
% L = spdiags(sum(Adj,2),0, N, N) - Adj;

% W-S random graphs
%N = 2^10;   % nodes
%K = 5;     % average degree (2*K)
%beta = 0;  % beta = 0 is a ring lattice, and beta = 1 is a random graph
%[~, L] = WattsStrogatz(N, K, beta); % N^2


% B-A random graphs
% N = 2000;
% A = ba_net('N', N);
% G =graph(A, 'OmitSelfLoops');
% L = laplacian(G);

%base = 7;
%L = make_poorly_conditioned(L, base);

%L = A_2D;
%L = S;

toc;

%---------------------
% get size
%---------------------
N = size(L,2);

%---------------------
% right hand side
%---------------------
%b = zeros(N,1);
%b = b_2D;
%exact_x = rand(N,1) + (1:N)';
%b = L*exact_x;
%b =zeros(N,1);
%b(2) = 1;
%b(3) = -1;
%b = buv;
%b = bN;
b = f(1:298431);

%b = b - sum(b)/N;

%---------------------
% AMG parameters
%---------------------
amgParam.print_level = 1; % how much information to print when using AMG solve only
                          % 0: print nothing | positive number print information  
% setup phase parameters
amgParam.max_level = 20; % maximal number of level in AMG
amgParam.coarsest_size = 100; % size of the coarest level

amgParam.strong_connection = 0.0; 
amgParam.agg_type = 'HEC';  % 'HEC': heavy edge coarsening | 'MIS': maximial indepent set coarsening
                                                    %| 'MWM': maximal weighted matching
amgParam.agg_radius = 1; % radius of the aggregates, only used when amgParam.agg_type = 'MIS'

% solve pahse parameters
amgParam.cycle_type = 'V'; % 'V': V-cycle | 'W': W-cycle | 'nV': n-fold V cycle | 'K': K-cycle (nonlinear AMLI-cycle)

amgParam.coarse_it = 2; % Number of iterations on coarse levels. Used when amgParam.cycle_type = 'nV' or 'K'
amgParam.coarse_krylov_type = 'GCG'; % Type of Krylove method when amgParam.cycle_type = 'K'
                                     % 'GCG': generalized CG | 'FGRMES': flexible GMRes 
                                     
amgParam.smooth_type = 'GS';   % 'Rd': Richardson | 'Jcb': Jacobi | 'GS': Gauss-Seidel | 'K1': Kaczmarz_AAt | 'K12': Kaczmarz_AtA
amgParam.n_presmooth = 1; % number of presmoothing      
amgParam.n_postsmooth =1; % number of postsmoothing
amgParam.smooth_omega = 1.2; % weight for Richardson smoother;

amgParam.prec_max_it = 1; % when AMG is used as preconditioner, how cycles will be applied
amgParam.max_it = 100;  % when AMG is used as standalone solver, maximal number of iterations that is allowed
amgParam.tol = 1e-8;    % when AMG is used as standalone solver, tolerance for the reletive residual

amgParam.ILU_level = 0;  % how many levels use ILU smoother (from fine to coarse), 0 means no ILU smoother
% ILU parameters (only used when amgParam.ILU_level > 0)
amgParam.ILU_type = 'crout';  % nofill, crout, ilutp
amgParam.droptol = 0.01; % drop tolerance for ILU when ILU_type = 'crout' or 'ilutp'
% do NOT touch the following three
amgParam.milu = 'off';  % row, col, off
amgParam.udiag = 1;
amgParam.thresh = 1;

%---------------------
% Iterative parameters
%---------------------
iterParam.solver_type = 'FGMRES'; % 'SL': GS as standalone solver | 'AMG': AMG as standalone solver | 'CG': conjugate gradiant | 'FGMRES': flexible GMRes
iterParam.prec_type = 'AMG'; % 'NULL': no preconditioner | 'AMG': AMG precondtioner
iterParam.print_level = 1; % how much information to print for iterative solver
                           % 0: print nothing | positive number print information 

iterParam.max_it = 100; % maximal number of iterations that is allowed
iterParam.tol = 1e-6;  % Tolerance for relative residual 
iterParam.restart = 100; % restart number for flexible GMRes

%---------------------
% setup phase
%---------------------
amgData = [];
if ( strcmp(iterParam.solver_type, 'AMG') || strcmp(iterParam.prec_type, 'AMG') )
    amgData = AMG_Setup(L, amgParam);
    %Lamg = L;
    %Lamg(end-23:end,end-23:end) =speye(24,24);
    %amgData = AMG_Setup(Lamg, amgParam);
end

%---------------------
% solve pahse
%---------------------
x = zeros(N,1); %initial guess 
%x = rand(N,1);

switch iterParam.solver_type
    
    case 'SL', % use simple linear iterative solver
        x = Simple_Solve(L, b, x, iterParam);

    case 'AMG', % use AMG directly
        x = AMG_Solve(amgData, b, x, amgParam);
        
    case 'LSQR' % use LSQR solver                                                                        
        % construct incidence matrix                                                                     
                                                                                                         
        L_ltrig = -tril(L,-1);                                                                           
        [s,t,weights] = find(L_ltrig);                                                                   
        G = graph(s,t,weights);                                                                          
        B = incidence(G);                                                                                
        m = length(G.Edges.Weight);                                                                      
        W_half = spdiags(G.Edges.Weight.^(1/2), 0, m, m);                                                
        %B_inv = pinv(full(B));                                                                          
        temp = rand(size(B,2),1);                                                                        
                                                                                                                                                                                   
        W_B = W_half * B';
        LSQR_start = tic;

        [sol,flag,relres,iter, resvec, lsvec] = lsqr(W_B, W_half * temp, 3e-12, 3000);
        residual = norm(L*sol - W_B'*W_half*temp)/norm(W_B'*W_half*temp)
        %[sol,flag,relres,iter, resvec, lsvec] = lsqr(B', omega, 1e-9, 3000);  % For MovieLens           
        %residual = norm(B*B'*sol - B*omega)/norm(B*omega)  % For MovieLens
        iter
        LSQR_duration = toc(LSQR_start) 
        
%    case 'LSRN' % use LSRN solver                                                                        
%         % construct incidence matrix                                                                     
%         L_ltrig = -tril(L,-1);                                                                           
%         [s,t,weights] = find(L_ltrig);                                                                   
%         G = graph(s,t,weights);                                                                          
%         B = incidence(G);                                                                                
%         m = length(G.Edges.Weight);                                                                      
%         W_half = spdiags(G.Edges.Weight.^(1/2), 0, m, m);                                                
%         temp = rand(size(B,2),1);                                                                        
%                                                                                                          
%         W_B = W_half * B';
%         x_init = rand(size(B,2),1);
%         [ sol, iter, duration, flag, relres,lsqr_res ] = LSRN( W_B, W_half * temp, 1e-11, 3000 );
%         %[ sol, iter, duration, flag, relres,lsqr_res ] = LSRN(B', omega, 1e-9, 3000, x_init);           
%        % lsqr_res                                                                                         
%         residual = norm(L*sol - W_B'*W_half*temp)/norm(W_B'*W_half*temp)                                 
%         %residual = norm(B*B'*sol - B*omega)/norm(B*omega)
%         duration
%         iter

    otherwise, % use preconditioned Krylov methods
        x = Krylov_Solve(L, b, x, iterParam, amgParam, amgData);
        
end

% pause;
% 
% iterParam.prec_type = 'NULL';
% x = zeros(N,1); %initial guess 
% x = Krylov_Solve(L, b, x, iterParam, amgParam, amgData);
