function [invAdiag, invA_diag] = invert_block_diag_matrix(A, block_size)

% get size
n = size(A,1);

% number of blocks
n_blocks = n/block_size;

% diagonal of inverse of A
invA_diag = cell(n_blocks,1);

% main loop
for i=1:n_blocks
   
    % index
    idx = ((i-1)*block_size+1):i*block_size;
    
    % get block
    block = A(idx, idx);
    
    % compute the inverse
    invA_diag{i} = inv(block);
    
end

invAdiag = blkdiag(invA_diag{:});

end

