function [ Q ] = generate_rand_proj( k, d, type )
% generate random projection
%
% INPUT:
%       k:      low dimension
%       d:      high dimension
%
% OUTPUT:
%       Q:      random projection matrix
%

switch type
    
    case 'gaussian' % random projection from normal distribution
        Q = 1/sqrt(k).*randn(k, d);
    
    case 'sparse'  % sparse random projection
        temp = rand(k,d);
        [row_p, col_p] = find(temp > 1-1/6);
        [row_n, col_n] = find(temp < 1/6);
        
        n_p = length(row_p);
        n_n = length(row_n);
        
        Q = sparse([row_p; row_n], [col_p; col_n], [sqrt(3/k)*ones(n_p,1); -sqrt(3/k)*ones(n_n,1)], k, d);
     
    otherwise  % dense random projection
        temp = rand(k,d);
        [row_p, col_p] = find(temp >= 0.5);
        [row_n, col_n] = find(temp < 0.5);
        
        n_p = length(row_p);
        n_n = length(row_n);
        
        Q = sparse([row_p; row_n], [col_p; col_n], [1/sqrt(k)*ones(n_p,1); -1/sqrt(k)*ones(n_n,1)], k, d);
        
end
        

end

