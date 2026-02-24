function [x] = kaczmarz_AAt(A, b, x, nsmooth, randomized)
% Kaczmarz algorithm for consistent system
% randomized takes in 1 and 0 to switch on/off the randomization
% @ Junyuan Lin & Xiaozhe Hu, Tufts University

%------------------------------------------
% get size
%------------------------------------------
[M, N] = size(A);

%------------------------------------------
% Calculate discrete CDF
%------------------------------------------
% tmp = A.^2;
% pm = sum(tmp, 2) / sum(tmp(:));
% F = cumsum(pm);

%------------------------------------------
% Calculate N random ordering
%------------------------------------------
% Choose row of A at random, with probability 
% proportional to norm(A(rp, :))^2
if randomized
    %rows = sum(repmat(F,1,N) < repmat(rand(1,N),M,1))+1;
    rows = randperm(N);
else
    rows = (1:N)';
end
a = A(rows,:);

% get lower triangular
DL = tril(a*a');

%------------------------------------------
% main loop
%------------------------------------------
for i = 1:nsmooth
        
        % iteration       
        %x = x + ( (b(rows) - a'*x) / (a' * a) )*a;
        x = x + a' * (DL\(b(rows) - a*x));

end

end


