function [ amgData ] = AMG_Setup_block( Af, amgParam )
% Setup phase for AMG method for matrices have block structure 
% (inseady each "entry" is a number, it is a nb by nb small matrix)
%
% @ Xiaozhe Hu, Tufts University

%----------------
% local variable
%----------------
print_level = amgParam.print_level;

max_level = amgParam.max_level;
coarsest_size = amgParam.coarsest_size;

agg_type = amgParam.agg_type;
agg_radius = amgParam.agg_radius;

nb = amgParam.nb;

level = 1;

%----------------
% AMG information
%----------------
P = cell(max_level,1);
R = cell(max_level,1);
A = cell(max_level,1);
N = cell(max_level,1);
DL = cell(max_level,1);
DU = cell(max_level,1);
D = cell(max_level,1);
Dinv = cell(max_level,1);

%----------------
% finest level
%----------------
A{1} = Af;
N{1} = size(Af,1);
[D{1}, DL{1}, DU{1}] = prepare_block_smoother(A{1}, nb);
[Dinv{1}] = invert_block_diag_matrix(A{1}, nb);

%----------------
% main loop
%----------------
if print_level >0 
    fprintf('----------------------------------------------------\n');
    fprintf('              Calling AMG setup    \n');
    fprintf('----------------------------------------------------\n');
end

setup_start = tic;

while ((level < max_level) && (N{level} > coarsest_size) )
        
    %----------------------------------------------------------------------
    % condense to a small matrix so coarsening algorithms can be performed
    %----------------------------------------------------------------------
    idx = 1:nb:N{level};
    As = A{level}(idx,idx);
    
    %----------------
    % strong connection
    %----------------
    As = construct_strong_connection( As, amgParam );
        
    %----------------
    % form aggregation
    %----------------
    switch agg_type
       
        case 'HEC'
            [ aggregation, num_agg ] = coarsening_agg_HEC(As);
        case 'MIS'
            [ aggregation, num_agg ] = coarsening_agg_MIS( As, agg_radius );
            %[ aggregation, num_agg ] = my_aggregation( As ); % Dylan
            %[ aggregation, num_agg ] = construct_aggregates(A{level}); % Duc
        case 'MWM'
            flag = 0;
            [ aggregation, num_agg ] = coarsening_agg_MWM( As, flag );
        otherwise
            disp('Wrong aggregation tyep!');
            disp('Use default type.')
            %[ aggregation, num_agg ] = coarsening_agg_MIS( As, agg_radius );
            [ aggregation, num_agg ] = coarsening_agg_HEC(As);
      
    end
    
    clear As;
    
    %----------------
    % generate prolongation
    %----------------
    [ Ps ] = generate_unsmoothed_P(aggregation, num_agg);
    P{level} = kron(Ps, speye(nb,nb));
        
    %----------------
    % generate restriction
    %----------------
    R{level} = P{level}';
    
    %----------------
    % compute coarse grid matrix
    %----------------
    A{level+1} = R{level}*A{level}*P{level};
    N{level+1} = size(A{level+1},1);
    
    %----------------
    % extra information
    %----------------    
    [D{level+1}, DL{level+1}, DU{level+1}] = prepare_block_smoother(A{level+1}, nb);
    [Dinv{level+1}] = invert_block_diag_matrix(A{level+1}, nb);
    
    %---------------
    % update
    %---------------
    if N{level+1}/N{level} >= 0.9
        amgParam.strong_connection = amgParam.strong_connection/8;
    elseif N{level+1}/N{level} >= 0.5 && N{level+1}/N{level} < 0.9
        amgParam.strong_connectionn = amgParam.strong_connection/4;
    end
    level = level+1;
        
end

setup_duration = toc(setup_start);

% make sure the coarsest level is big enough for eigenvalue problems
if (N{level}<=amgParam.number_eigen)
    level = level - 1;
end

% construct the data structure
amgData = struct('P',P,'R',R,'A',A,'DL',DL,'DU',DU,'D',D,'Dinv',Dinv,...
    'N',N,'max_level',level,'number_eigen',amgParam.number_eigen,'nb',nb);

amgData = amgData(1:level);

% print information
if print_level > 0
   
    total_N = 0;
    total_NNZ = 0;
    
    fprintf('----------------------------------------------------\n');
    fprintf(' # Level |   # Row   |   # Nonzero  | Avg. NNZ/Row |\n');
    fprintf('----------------------------------------------------\n');
    
    for i = 1:level
            
            total_N = total_N + N{i};
            total_NNZ = total_NNZ + nnz(A{i});
            
            fprintf('   %2d    |%9d  | %10d   |   %7.3f    |\n', i, N{i}, nnz(A{i}), nnz(A{i})/N{i});
    end
    
    fprintf('----------------------------------------------------\n');
    
    fprintf(' Grid complexity: %0.3f | Operator complexity: %0.3f \n', total_N/N{1}, total_NNZ/nnz(A{1}));
    
    fprintf('----------------------------------------------------\n');
    
end

% print cputime
if print_level > -1

    fprintf('----------------------------------------------------\n');
    fprintf('        AMG setup costs %f seconds\n', setup_duration);
    fprintf('----------------------------------------------------\n');

end

end

