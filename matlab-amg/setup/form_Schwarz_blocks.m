function [schwarzData] = form_Schwarz_blocks(seeds, A)
%

% get size
N_blocks = length(seeds);

% get cell flag
flag_cell = iscell(seeds);

% initialize 
blocks = cell(N_blocks,1);

% form blocks
for i = 1:N_blocks
    
    % nodes on the current radius
    if flag_cell
        current_seeds = seeds{i};
    else
        current_seeds = seeds(i);
    end
    
    % neighboring DoFs 
    [~,jj,~] = find(A(current_seeds,:));
    
    % get the block
    blocks{i} = unique(jj);
    
end

% block smoother data
schwarzData = struct('N_blocks', N_blocks, 'blocks', blocks);

end

