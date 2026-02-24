function [L_minus, L_M_rel, L_tilde, ind_reorder,P_mat,layer_ind] = build_MSP(L_1, alpha, level,agg_type)
% build multilevel subgraph preconditioner

% Inputs: L_1 - graph laplacian for the original graph
%         alpha - scalar to control the condition number of L_tilde*L_M_rel (denoted as mu in the paper) 
%         level - # of levels of the preconditioner
%         agg_type - aggregation type

% Outputs: L_T - graph laplacian of the MSP
%          L_M_rel - graph laplacian of the positive subgraph
%          L_tilde - graph laplacian of the expanded graph from L_1
%          P_mat - prolongation matrix for generating b_tilde
%          layer_ind - layer indices for rearranging the matrices
%clear;
%clc;


n = size(L_1,1);
e = ones(n,1);

switch agg_type
    
    case 'matching'
        temp_P = cell(level, 1);
        temp_P{1,1} = spdiags(e,0,n,n);
        temp_level = 2;
        num_agg_total = cell(level+1,1);
        num_agg_total{1} = 0;
        num_agg_total{2} = n;
        %nextLevelP
        while temp_level < level+1 
            if mod(temp_level,2)==0
                nextLevelP = sparse(size(L_1,1)/(2^(temp_level-2)), size(L_1,1)/(2^(temp_level-1)));
                for i = 1:size(nextLevelP,2)
                    nextLevelP(2*i-1,i)=1;
                    nextLevelP(2*i,i) =1;
                end
            else
                nextLevelP = sparse(size(L_1,1)/(2^(temp_level-2)), size(L_1,1)/(2^(temp_level-1)));
                for i = 1:size(nextLevelP,2)
                    nextLevelP(2*i-1,i)=1;
                    nextLevelP(2*i,i) =1;
                end
        %   % for grids, make sure the aggregations are the nodes on each square grid
        %         num_parts = (level-temp_level-1)/2;
        %         num_parts = 2^(num_parts+1);
        %         nextLevelP = [];
        %         nextLevelP_temp = sparse(size(L_1,1)/(2^(temp_level-2))/num_parts, size(L_1,1)/(2^(temp_level-1))/num_parts);
        %         for i = 1:size(nextLevelP_temp,2)
        %                 nextLevelP_temp(i,i)=1;
        %                 nextLevelP_temp(i+num_parts,i)=1;
        %         end
        %         for k = 1:num_parts
        %             nextLevelP = blkdiag(nextLevelP,nextLevelP_temp);
        %         end
            end
            temp_P{temp_level, 1} = nextLevelP;
            num_agg_total{temp_level+1} = num_agg_total{temp_level} + size(nextLevelP,2);
            temp_level = temp_level + 1;
        end
       
    case 'MWM'
        As = L_1;
        l = 2;
        temp_P{1,1} = spdiags(e,0,n,n);
        num_agg_total = cell(level+1,1);
        num_agg_total{1} = 0;
        num_agg_total{2} = n;
        while l < level+1
            [ aggregation, num_agg ] = coarsening_agg_MWM( As, 1 );
            [ temp_P{l,1} ] = generate_unsmoothed_P(aggregation, num_agg);
            As = temp_P{l,1}'*As*temp_P{l,1};
            num_agg_total{l+1} = num_agg_total{l} + num_agg;
            l = l+1;
        end
end
%P stores all the prolongation operators P_1^2, P_1^3 etc. in each cell
P = cell(level, 1);
for i = 1:level
    temp_P_mat = 1;
    for j = 1:i
        temp_P_mat = temp_P_mat*temp_P{j,1};
    end
    P{j, 1} = temp_P_mat;
end

%%hard-coded for checking if three layer works - can delete
% block_L = cell(level, 1);
% for i = 1:level
%     block_L{i, 1} = P{i, 1}'*L_1*P{i, 1};
% end
% L_1_tilde1 = blkdiag(block_L{1, 1},block_L{2, 1},block_L{3, 1}); 


%build L_tilde to take the block diagonals of L_G
n_tilde = num_agg_total{end,1};
temp_mat = sparse(n_tilde,n);
%temp_mat(1:n,:) = P{1, 1}';
%L_1_tilde = temp_mat*L_1*temp_mat';
L_1_tilde = sparse(n_tilde,n_tilde);
for i = 1:level
    temp_mat = sparse(n_tilde,n);
    starting = num_agg_total{i}+1;
    temp_mat(starting:num_agg_total{i+1}, :) = (-1)^(i-1)*P{i, 1}';
    L_1_tilde = L_1_tilde + temp_mat*L_1*temp_mat';
end



%hard-coded for checking
% for i = 1:size(P_1_2,2)
%     P_1_2(2*i-1,i)=1;
%     P_1_2(2*i,i)=1;
% end
% P_2_3 = zeros(size(L_1,1)/2, size(L_1,1)/4);
% for i = 1:size(P_2_3,2)
%     P_2_3(2*i-1,i)=1;
%     P_2_3(2*i,i)=1;
% end
% P_1_3 = P_1_2*P_2_3;
% unstructured P
%P_1 = [1 0 0 0; 0 1 1 1]';

% P_2 = zeros(size(P_1,2), size(P_1,2)/2);
% for i = 1:size(P_2,2)
%     P_2(2*i-1,i)=1;
%     P_2(2*i,i)=1;
% end
% 
% 
% P_2 = P_1*P_2;
%Id = eye(size(L_1,1));


%L_tilde = [-P{1,1}';-P{2,1}'; -P{3,1}']*L_1*[-P{1,1}',-P{2,1},-P{3,1}];

P_mat = sparse(n_tilde,n);
%P_mat(1:n,:) = P{1, 1}';
for i = 1:level
    starting = num_agg_total{i}+1;
    P_mat(starting:num_agg_total{i+1}, :) = (-1)^(i-1)*P{i, 1}';
end
L_tilde = P_mat*L_1*P_mat';


%%%%%%%%%%%%%%%%%%%%%%
%   open if needed
% check diameter for L_M_rel_unscaled
%L_M_rel_unscaled = L_tilde;
%pos = find(L_M_rel_unscaled>0);
%L_M_rel_unscaled(pos) = 0;
%diagonals = sum(L_M_rel_unscaled,2); 
%L_M_rel_unscaled = L_M_rel_unscaled - diag(diagonals);
%L_M_rel_unscaled = (L_M_rel_unscaled + L_M_rel_unscaled')/2;

%check L_M_rel_unscaled's diameter
%diameter_unscaled=diameter(L_M_rel_unscaled)
%%%%%%%%%%%%%%%%%%%%%%%%%%

L_minus = L_tilde;
L_minus = L_minus - spdiags(diag(L_tilde), 0, size(L_tilde, 1), size(L_tilde, 1));
pos = find(L_minus<0);
L_minus(pos) = 0;
diagonals = sum(L_minus);
L_minus = -L_minus + diag(diagonals);
%r = qr(L_1_tilde,0); %for finding pinv on sparse L_1_tilde
%lambda_max = max(eigs(L_minus*(r\(r'\L_1_tilde')))) %pinv of L_1_tilde=r\(r'\A')
if strcmp(alpha,"max")
    
    %estimate lambda(max(eig(pinv(full(L_1_tilde))*L_minus))) by
    %lambda_max(L_minus)*lambda_min(L_1_tilde)
%     n=length(L_minus); % Size of initial eigenvector.
%     u1=zeros(length(L_minus),1); % The initial choice of eigenvector.
%     u1(1) = 1;
%     u2 = u1;
%     for iter_init = 1:10  %Calculating the greatest eigenvalue and the corresponding eigenvector.
%        v1=L_minus*u1; 
%        u1=v1/max(abs(v1));
%        v2 = L_1_tilde\u2;
%        u2 = v2/norm(v2,2);
%     end
%     lambda_max = max(abs(v1))*min(abs(v2));
    [v_L_1_tilde,eigenvals_L_1_tilde] = eig(full(L_1_tilde));
    Z = sparseinv(L_1_tilde);
    eigenvals_Z = eigs(Z,5,'smallestabs');
    [v_L_minus,eigenvals_L_minus] = eig(full(L_minus));
    lambda_max = max(eigs(L_minus*Z));
%     lambda_max = max(eig(pinv(full(L_1_tilde))*L_minus))
%     alpha = 1/(2*lambda_max)
end
%
%build steiner graph L_tilde
%n=2;
%alpha = 1/2
%alpha = 1;
%Id = eye(size(L_1,1));
P_mat = sparse(n_tilde,n);
%P_mat(1:n,:) = P{1, 1}';
for i = 1:level
    starting = num_agg_total{i}+1;
    P_mat(starting:num_agg_total{i+1}, :) = (-1*alpha)^(i-1)*P{i, 1}';
end
L_tilde = P_mat*L_1*P_mat';
%L_tilde = [Id;-alpha*P_1_2'; -alpha*P_1_3']*L_1*[Id,-alpha*P_1_2,-alpha*P_1_3];
%L_tilde = [Id;-alpha*P_1_2']*L_1*[Id,-alpha*P_1_2];
L_tilde = (L_tilde+L_tilde')/2;

%%%%%%%%%%%%%%%%%%%%%%
%   open if needed
%%L_tilde = [3 -4 1; -4 5 -1; 1 -1 0];
% eig_L_tilde=eig(L_tilde);
% eig_L_tilde_pos = find(eig_L_tilde>1e-5);
% max_eig_L_tilde = max(eig_L_tilde(eig_L_tilde_pos));
% min_eig_L_tilde = min(eig_L_tilde(eig_L_tilde_pos));
% cond_L_tilde = max_eig_L_tilde/min_eig_L_tilde
%%%%%%%%%%%%%%%%%%%%%%

%L_tilde = L+ - L-
% %construct M relative matrix of L_tilde (L_tilde+)
L_M_rel = L_tilde;

% for i = size(L_1,1) + 1 : size(L_tilde,1)
%     for j = 1:size(L_1,1)
%         if L_M_rel(i,j)>0
%             L_M_rel(i,i) = L_M_rel(i,i) + L_M_rel(i,j);
%             L_M_rel(j,j) = L_M_rel(j,j) + L_M_rel(i,j);
%             L_M_rel(i,j) = 0;
%             L_M_rel(j,i) = 0;
%         end
%     end
% end

%L_M_rel = L_M_rel - spdiags(diag(L_tilde), 0, size(L_tilde, 1), size(L_tilde, 1));
pos = find(L_M_rel>0);
L_M_rel(pos) = 0;
%c = 1; %c is a scaling, c>1
%L_M_rel(1,2) = -c;
%L_M_rel(2,1) = -c;
diagonals = sum(L_M_rel,2); 
L_M_rel = L_M_rel - diag(diagonals);
L_M_rel = (L_M_rel + L_M_rel')/2;

%check L_M_rel's diameter
%diameter=diameter(L_M_rel)

%eig_H = eig(full(L_M_rel));

%L_M_rel = c*L_M_rel;
% eig_L_M_rel = eig(L_M_rel);
% eig_L_M_rel_pos = find(eig_L_M_rel>1e-5);
% max_eig_L_M_rel = max(eig_L_M_rel(eig_L_M_rel_pos));
% min_eig_L_M_rel = min(eig_L_M_rel(eig_L_M_rel_pos));
% cond_L_M_rel = max_eig_L_M_rel/min_eig_L_M_rel

% %construct L_tilde^- matrix of L_tilde (Laplacian of negatively weighted edges)
L_minus = L_tilde;
L_minus = L_minus - spdiags(diag(L_tilde), 0, size(L_tilde, 1), size(L_tilde, 1));
pos = find(L_minus<0);
L_minus(pos) = 0;
diagonals = sum(L_minus);
L_minus = -L_minus + diag(diagonals);
%eig_L_minus = eig(L_minus)


%% build tree for L_M_rel
% L_T = L_M_rel;
% %since the aggregations are horizontal matches, take the vertical lines on
% %L_1 out
% for i = 1: (size(L_1,1)-n)
%     L_T(i, i+n) = 0;
%     L_T(i+n, i) = 0;
% end
% %take out horizontal edges that are out of the aggregations
% for i = 1: size(L_1,1)/2
%     L_T(2*i, 2*i+1) = 0;
%     L_T(2*i+1, 2*i) = 0;
% end
% %take out some edges between 1 & 2 layers
% rep = [];
% for i = 1:n/2
%     rep = [rep, 2*i-1:2*n:size(L_1,1)];
% end
% rep = [rep, rep+n+1];
% n_l2 = (size(L_1,1)+1): size(L_1,1)+size(L_1,1)/2;
% for i = 1:length(rep)
%     L_T(rep(i), n_l2) = zeros(length(n_l2),1);
%     L_T(n_l2, rep(i)) = zeros(length(n_l2),1);
% end

%take out MWST from the second layer
% second_layer = L_T(n_l2, n_l2);
% [ second_layer_L_T, Cost ] =  UndirectedMaximumSpanningTree (second_layer);
% L_T(n_l2, n_l2) = second_layer.*second_layer_L_T;
L_T = sparse(size(L_M_rel,1), size(L_M_rel,2));
layer_ind = cell(level,1);
%layer_ind{1,1} = 1:num_agg_total{1};
for i = 1:level
    starting = num_agg_total{i}+1;
    layer_ind{i,1} = starting:num_agg_total{i+1};
end
%keep only the edges in the matchings on each layer (upper triangular)
ind_reorder = [];
for i = 2:level
    agg_edge = temp_P{i,1};
    for j = 1:size(agg_edge,2)
        ind = find(agg_edge(:,j))+num_agg_total{i-1}; % since it is not exactly a matching, so we need to add a path made of nodes in the aggregations
        L_T(ind, ind) = L_M_rel(ind, ind);
        %for k = 1:length(ind)-1
        %    L_T(ind(k), ind(k+1)) = L_M_rel(ind(k), ind(k+1));
        %end
        ind_reorder = [ind_reorder; ind];
    end
end
ind_reorder = [ind_reorder; layer_ind{level,1}'];
%ind_reorder = [];
% for i = 2:level
%     layer_mat = (ones(length(layer_ind{i,1}),1)*layer_ind{i-1,1})';
%     layer_mat = layer_mat.*temp_P{i,1};
%     for j = 1:length(layer_ind{i,1})
%         [~,~,ind] = find(layer_mat(:,j)); % should be a vector of length 2 since it's a matching
%         L_T(ind(1), ind(2)) = L_M_rel(ind(1), ind(2)) ;
%         %ind_reorder = [ind_reorder; ind];
%     end
% end
%ind_reorder = [ind_reorder; layer_ind{level,1}'];
% copy edges between two consecutive layers from L_M_rel
for i = 1:level-1
    L_T(layer_ind{i,1},layer_ind{i+1,1}) = L_M_rel(layer_ind{i,1},layer_ind{i+1,1});
end
L_T = L_T + L_T';
% copy the top level graph Laplacian in case if the top level is not just
% one aggregation
L_T(layer_ind{level,1}, layer_ind{level,1}) = L_M_rel(layer_ind{level,1}, layer_ind{level,1});


% take out some edges between first and second layers
% rep = [];
% for layer_count = 1:level-1
%     lines = floor(sqrt(length(layer_ind{i,1})));
%     for i = 1:lines/2
%         rep = [rep, 2*i-1:2*lines:length(layer_ind{i,1})];
%     end
%     rep = [rep, rep+lines+1];
% end

%% spectral analysis
%preconditioner
% ST = diag(diag(L_1)); %ST preconditioner here we need to only take the diagonals
% up_st = ST*(-1)*alpha*P_1_2;
% down_st = (-1)*alpha*P_1_2'*ST;
% L_T = [zeros(size(L_1)) up_st;down_st (-1)*alpha*P_1_2'*(L_1)*(-1)*alpha*P_1_2];
% %make L_T a graph Laplacian
% L_T = spdiags(zeros(size(L_T,1),1), 0, L_T);
% L_T = spdiags(-sum(L_T)', 0, L_T);
% % %L_T = [ST (-1)*L_1*P_1;(-1)*P_1'*L_1 (-1)*P_1'*(L_1+ST)*(-1)*P_1];
% % %L_T = [ST up_st;down_st (-1)*P_1'*L_1*(-1)*P_1];
% eig_L_T = eig(pinv(full(L_T))*L_tilde);
% eig_L_T_pos = find(eig_L_T>1e-5);
% max_eig_L_T = max(eig_L_T(eig_L_T_pos));
% min_eig_L_T = min(eig_L_T(eig_L_T_pos));
% % cond_L_T = max_eig_L_T/min_eig_L_T

% lamda(L+/L-)
% [v_L_pm,eig_L_pm] = eig(pinv(full(L_M_rel))*L_minus);
% eig_L_pm_pos = find(eig_L_pm>1e-5);
% max_eig_L_pm = max(eig_L_pm(eig_L_pm_pos));
% min_eig_L_pm = min(eig_L_pm(eig_L_pm_pos));
% cond_L_pm = max_eig_L_pm/min_eig_L_pm;

%%%%%%%%%%%%%%%%%%%%%%
%   open if needed
%lamda(L_tilde/L_M_rel)
% eig_L_M_rel = eig(pinv(full(L_M_rel))*L_tilde);
% eig_L_M_rel_pos = find(eig_L_M_rel>1e-5);
% max_eig_L_M_rel = max(eig_L_M_rel(eig_L_M_rel_pos));
% min_eig_L_M_rel = min(eig_L_M_rel(eig_L_M_rel_pos));
% one_over_min = 1/min_eig_L_M_rel;
% cond_L_M_rel = max_eig_L_M_rel/min_eig_L_M_rel;
%est_cond = lambda_max*alpha^(2-2*level)


% lamda(L_M_rel/L_T)
%make L_T a graph Laplacian
L_T_n = size(L_T,1);
d = diag(L_T);
L_T = (L_T-spdiags(d,0,L_T_n,L_T_n));
%L_T = -spones(L_T);
d = sum(L_T);
if size(d,1)==1
   d = d';
end 
L_T = L_T-spdiags(d,0,L_T_n,L_T_n);


%[v_Lt, eig_Lt] = eig(full(L_T));
% make L_T orthogonal to its nullspace
%L_T_0 = 
%[v_Lm, eig_Lm] = eig(full(L_M_rel));
% make L_M_rel orthogonal to its nullspace

%[v_test, eig_L_T] = eig(pinv(full(L_T))*L_M_rel);
%eig_L_T = eig(pinv(full(L_T))*L_M_rel);
%[v,eig_val] = eigs(L_M_rel,L_T,509,'sm');
%eig_L_T_pos = find(eig_L_T>1e-5);
% max_eig_L_T = max(eig_L_T(eig_L_T_pos));
% min_eig_L_T = min(eig_L_T(eig_L_T_pos));
% %one_over_min = 1/min_eig_L_T
% cond_L_T = max_eig_L_T/min_eig_L_T

%%%%%%%%%%%%%%%%%%%%%%
%   open if needed
%total stretch
%[total_stretch] = calculate_stretch(L_M_rel,L_T)


%multigrid
% L_tg = [ST zeros(size(L_1,1), size(P_1, 2));-P_1'*L_1 P_1'*L_1*P_1];
% eig_L_tg = eig(pinv(full(L_tg))*L_tilde);
% eig_L_tg_pos = find(eig_L_tg>1e-5);
% max_eig_L_tg = max(eig_L_tg(eig_L_tg_pos))
% min_eig_L_tg = min(eig_L_tg(eig_L_tg_pos))
% cond_L_tg = max_eig_L_tg/min_eig_L_tg


% 
% ST_2 = diag(diag(L_T));
% up_st = ST_2*(-1)*P_2;
% down_st = (-1)*P_2'*ST_2;
% L_T = [ST up_st;down_st (-1)*P_1'*ST*(-1)*P_1];


% % % %plot for original graph
% %figure
% L_tilde = (L_tilde + L_tilde')/2;
% G = graph(L_tilde,'omitself');
% %m_G = size(G.Edges,1)
% %plot(G,'Layout','force')
% plot(G,'EdgeLabel',-G.Edges.Weight);
% % % % % 
% % % % %plot for M relative graph
% figure;
% L_M_rel = (L_M_rel+L_M_rel')/2;
% G_M_rel = graph(L_M_rel,'omitself');
% m = size(G_M_rel.Edges,1);
% plot(G_M_rel,'Layout','force')
% %plot(G_M_rel,'EdgeLabel',-G_M_rel.Edges.Weight);
% % % 
% 
% figure;
% L_minus = (L_minus+L_minus')/2;
% G_minus = graph(L_minus,'omitself');
% plot(G_minus,'EdgeLabel',-G_minus.Edges.Weight);
% % 
%%plot for the preconditioner
% L_T = (L_T+L_T')/2;
% G_T = graph(L_T,'omitself');
% figure;
% plot(G_T, 'Layout','force')
%plot(G_T,'EdgeLabel',-G_T.Edges.Weight)
% pause
% calculating exact eigenvalues
% d_L_tilde = eig(L_tilde);
% d_L_M_rel = eig(L_M_rel)


%for calculating the upperbound of the preconditioner
% [ D_G,B_G ] = construct_BD( L_tilde );
% [ D_T,B_T ] = construct_BD( L_T );
% %W = pinv(B_T)*B_G;
% [ W ] = construct_W( B_G, B_T);
% 
% diag_arr = diag(D_T);
% diag_half = diag_arr.^(-1/2);
% 
% % D_T^(-1/2)*W*D_G*W'*D_T^(-1/2)
% W_tilde = diag(diag_half)*W*D_G*W'*diag(diag_half);
% 
% %exact
% eig_W_tilde=eig(W_tilde);
% eig_W_tilde_pos = find(eig_W_tilde>1e-5);


% Upperbound
%norm_2 = norm(W_tilde, 2);
%norm_2_sq = norm_2^2
%max_W_t = max(eig_W_tilde(eig_W_tilde_pos))
% max_e_sq = max(eig_W_tilde(eig_W_tilde_pos))^2
% norm_F = norm(W_tilde, 'fro');
% norm_F_sq = norm_F^2
% norm_1i_est = norm(W_tilde, 1) * norm(W_tilde, inf)
% upbound_lim_pf = (1 + sqrt(max(deg)/min(deg)))^2

%lowerbound
%min_W_t = min(eig_W_tilde(eig_W_tilde_pos))
% lowerbound_lim_pf =  1/(1 + 4*max(sum(P_1)))
% 
% %bound
%exact_bd = max(eig_W_tilde(eig_W_tilde_pos))/min(eig_W_tilde(eig_W_tilde_pos))
% limbd_pf = upbound_lim_pf/lowerbound_lim_pf 





% max_eval_approx = max(max(D))
% for i = 1:size(W_tilde,1)
%     norm_arr(i) = norm(W_tilde(i,:),1);
% end
% [norm_max,ind] = max(norm_arr)
% upper_bound = norm(W_tilde,Inf)
% est_bound1 = 3*max(deg);
% est_bound2 = 2*size(L_1,1)+4;

% below is probably not right, but not sure where
% [ S ] = construct_S( B_G,B_T );
% diag_arr = diag(D_G);
% diag_arr = diag_arr.^(-1);
% mat2 = diag(diag_arr)*S*D_T*S';
% [~,D]=eig(mat2);
% temp_max = max(max(abs(D)));
% min_eval_approx = 1/temp_max


%for calculating the lowerbound of the preconditioner
% [ L_H,L_nng ] = construct_LH( L_tilde,L_T,size(L_1,1) );
% [ D_H,B_H ] = construct_BD( L_H );
% [ D_nng,B_nng ] = construct_BD( L_nng );
% [ M ] = construct_M( B_nng,B_H);
% 
% diag_arr = diag(D_H);
% diag_arr = diag_arr.^(-1);
% mat1 = diag(diag_arr)*M*D_nng*M';
% [~,D_min]=eig(mat1);
% min_eval_approx = max(max(D_min))
%min_eval_approx = 1/min_eval_approx




end

