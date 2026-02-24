function [x] = kaczmarz_AtA(A, b, x, nsmooth, randomized)
% Kaczmarz algorithm for consistent system
% randomized takes in 1 and 0 to switch on/off the randomization
% @ Junyuan Lin & Xiaozhe Hu, Tufts University

%------------------------------------------
% get size
%------------------------------------------
[~, N] = size(A);

%------------------------------------------
% Calculate discrete CDF
%------------------------------------------
% tmp = A'.^2;
% pm = sum(tmp, 2) / sum(tmp(:));
% F = cumsum(pm);

%------------------------------------------
% Calculate N random ordering
%------------------------------------------
% Choose row of A at random, with probability 
% proportional to norm(A(rp, :))^2
if randomized
    %cols = sum(repmat(F,1,N) < repmat(rand(1,N),N,1))+1;
    %cols = randi([1,N],N,1);
    cols = randperm(N);
else
    cols = (1:N)';
end

% get lower triangular
DL = tril(A(:,cols)'*A(:,cols));

%------------------------------------------
% main loop
%------------------------------------------
for i = 1:nsmooth
        
        % iteration       
        x = x +  DL\(A'*(b - A*x));

end

end


