function [D, DL, DU] = prepare_block_smoother(A, nb)
% get block diagonal, lower/upper triangular of a matrix with "entry-wise"
% block sturctureb
%
% @ Xiaozhe Hu, Tufts University 

% get size
n = size(A,1);

% allocate
D = sparse(n,n);
DL = sparse(n,n);
DU = sparse(n,n);

% size of each big block
blk_size = n/nb;

% main loop
for ib=1:nb
    
    % get row index
    idx_row = ib:nb:n;
    
    for jb=1:nb
       
        % get colume index
        idx_col = jb:nb:n;
        
        % get subblock
        subA = A(idx_row, idx_col);
        
        % get diagonal block
        subD = spdiags(diag(subA),0, blk_size, blk_size);
        [i,j,val] = find(subD);
        D = D + sparse((i-1)*nb+ib,(j-1)*nb+jb, val, n, n);
        
        % get lower triangular block
        subDL = tril(subA); 
        [i,j,val] = find(subDL);
        DL = DL + sparse((i-1)*nb+ib,(j-1)*nb+jb, val, n, n);

        % get upper triangular block
        subDU = triu(subA); 
        [i,j,val] = find(subDU);
        DU = DU + sparse((i-1)*nb+ib,(j-1)*nb+jb, val, n, n);
        
    end 
end

end