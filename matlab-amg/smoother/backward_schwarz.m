function [x] = backward_schwarz(A, b, x, nsmooth, schwarzData)
%

% get number of blocks
N_blocks = schwarzData.N_blocks;

% outter loop
for k = 1:nsmooth
    
    for i = N_blocks:-1:1
        
        % get the block
        current_block = schwarzData(i).blocks;
        
        % compute residual (needs a better way to do this -- XH)
        r = b - A*x;
        
        % solve locally 
        e_local = A(current_block,current_block)\r(current_block);
        
        % update
        x(current_block) = x(current_block) + e_local;
        
    end

end

