function [L] = make_poorly_conditioned(L, base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes in a graph Laplacian and makes it poorly conditioned
%by randomly altering the edge weights to be 10^b or 10^(-b) for some
%integer b. This is effectively just scaling the matrix from the left and
%right. 
%Input:
%   L - 
%       a graph Laplacian 
%Returns:
%   L -
%       a sparse matrix with the altered edge weights of L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% base = 12;
% large_num = 10^(base);
% small_num = 10^(-base);


%L will be symmetric so we will process it as a striculy upper triangular matrix
U = triu(L,1);

[i,j,edge_weights] = find(U);
edge_count = length(edge_weights);

% for k = 1:edge_count
%     
%     coin_flip = rand;
%     if coin_flip > .5
%         edge_weights(k) = large_num;
%     else
%         edge_weights(k) = small_num;
%     end
% end

% coin_flip = rand(edge_count, 1);
% flag = (coin_flip > .5);
% edge_weights(flag) = large_num;
% edge_weights(~flag) = small_num;

edge_weights = exp(base*(rand(edge_count, 1) - 0.5));

%remake L 

%make symmetric
new_L = sparse([i,j],[j,i],repmat(edge_weights,1,2));
degrees = sum(new_L); %compute row sums
L = diag(degrees) - new_L;


end
