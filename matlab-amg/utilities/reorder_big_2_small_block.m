function [Asmall, idx] = reorder_big_2_small_block(Abig, nb)
% reorder a matrix with big block structure to "entry-wise" block structure 
%
% @ Xiaozhe Hu, Tufts University 

% get size
n = size(Abig,1);

% size of each big block
n_block = n/nb;

% get idx
idx = (reshape((1:n), n_block, 3))';
idx = idx(:);

% reorder
Asmall = Abig(idx,idx);

end