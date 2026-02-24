function Abig = reorder_small_2_big_block(Asmall, nb)
% reorder a matrix with "entry-wise" block structure to big block structure
%
% @ Xiaozhe Hu, Tufts University 

% get size
n = size(Asmall,1);

% size of each big block
n_block = n/nb;

% get idx
idx = (reshape((1:n), 3, n_block))';
idx = idx(:);

% reorder
Abig = Asmall(idx,idx);

end

