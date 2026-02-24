function [kappa, T] = condition_estimate_CG(alpha,beta)
% using the CG coefficients to estimate the condition number
%
% @ Xiaozhe Hu, Tufts University

% get size
N = length(alpha);

% beta/alpha
temp = [0; beta./alpha(1:N-1)];

% diagonal entries
d0 = 1./alpha + temp;

% offdiaognal entries
d1 = sqrt(beta)./alpha(1:N-1);

% form T
T = spdiags([d1;0], -1, N, N) + spdiags(d0, 0, N, N) + spdiags([0;d1], 1, N, N);

% compute extreme eigenvalues of T
lambda_max = eigs(T,1,'lm');
lambda_min = eigs(T,1,'sm');

% compute the condition number
kappa = lambda_max/lambda_min;

end

