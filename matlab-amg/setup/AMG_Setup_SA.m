function [ amgData ] = AMG_Setup_SA( Af, amgParam )
% Setup phase for smoothed aggregation AMG (SAAMG) method
%
% @ Xiaozhe Hu, Tufts University

%----------------
% local variable
%----------------
print_level = amgParam.print_level;

max_level = amgParam.max_level;
coarsest_size = amgParam.coarsest_size;

Schwarz_level = amgParam.Schwarz_level;

ILU_level = amgParam.ILU_level;
if ILU_level > 0 
    ILU_setup.type = amgParam.ILU_type;  % nofill, crout, ilutp
    ILU_setup.droptol = amgParam.droptol;
    ILU_setup.milu = amgParam.milu;  % row, col, off
    ILU_setup.udiag = amgParam.udiag;
    ILU_setup.thresh = amgParam.thresh;
end


agg_type = amgParam.agg_type;
agg_radius = amgParam.agg_radius;

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
schwarzData = cell(max_level,1);
IL = cell(max_level, 1);
IU = cell(max_level, 1);
infNorm = cell(max_level, 1);

%----------------
% finest level
%----------------
A{1} = Af;
N{1} = size(Af,1);
I = speye(N{1},N{1});
DL{1} = tril(Af); % + 1e-8*I;
DU{1} = triu(Af); % + 1e-8*I;
D{1} = diag(Af); % + 1e-8;
Dinv{1} = 1./D{1};
infNorm{1} = norm(Af, 'inf');

%----------------
% main loop
%----------------
if print_level >0 
    fprintf('----------------------------------------------------\n');
    fprintf('              Calling SA AMG setup    \n');
    fprintf('----------------------------------------------------\n');
end

setup_start = tic;

while ((level < max_level) && (N{level} > coarsest_size) )
    
    %----------------
    % strong connection
    %----------------
    As = construct_strong_connection( A{level}, amgParam );
        
    %----------------
    % form aggregation
    %----------------
    switch agg_type
       
        case 'HEC'
            %[ aggregation, num_agg ] = coarsening_agg_HEC(A{level});
            [ aggregation, num_agg ] = coarsening_agg_HEC(As);
        case 'MIS'
            %[ aggregation, num_agg ] = coarsening_agg_MIS( A{level}, agg_radius );
            [ aggregation, num_agg ] = coarsening_agg_MIS( As, agg_radius );
            %[ aggregation, num_agg ] = my_aggregation( As ); % Dylan
            %[ aggregation, num_agg ] = construct_aggregates(A{level}); % Duc
        case 'MWM'
            flag = 1;
            [ aggregation, num_agg ] = coarsening_agg_MWM( As, flag );
            %[ aggregation, num_agg ] = coarsening_agg_MWM( A{level}, flag );
        otherwise
            disp('Wrong aggregation type!');
            disp('Use default type: heavy edge coarsening');
            %[ aggregation, num_agg ] = coarsening_agg_MIS( A{level}, agg_radius );
            %[ aggregation, num_agg ] = coarsening_agg_MIS( As, agg_radius );
            [ aggregation, num_agg ] = coarsening_agg_HEC(As);
      
    end
        
    %----------------
    % generate the tentative prolongation
    %----------------
    [ P{level} ] = generate_unsmoothed_P(aggregation, num_agg);
    
    %----------------
    % smooth the prolongation
    %----------------
    [ P{level} ] = generate_smoothed_P_Jacobi( P{level}, As );
        
    clear As;
        
    %----------------
    % generate the restriction
    %----------------
    R{level} = P{level}';
    
    %----------------
    % compute the coarse grid matrix
    %----------------
    A{level+1} = R{level}*A{level}*P{level};
    N{level+1} = size(A{level+1},1);
    
    %----------------
    % extra information
    %----------------
    if level <= Schwarz_level
        if (level == 1) && (isfield(amgParam, "Schwarz_blocks"))
            schwarzData{level} = struct('N_blocks', length(amgParam.Schwarz_blocks), 'blocks', amgParam.Schwarz_blocks);
        else
            seeds = 1:N{level};
            schwarzData{level} = form_Schwarz_blocks(seeds, A{level});
        end
    end 

    if level <= ILU_level
        [IL{level}, IU{level}] = ilu(sparse(A{level}), ILU_setup); 
    end
        
    I = speye(N{level+1},N{level+1});
    DL{level+1} = tril(A{level+1}); % + 1e-8*I;
    DU{level+1} = triu(A{level+1}); % + 1e-8*I;
    D{level+1} = diag(A{level+1}); % + 1e-8;
    Dinv{level+1} = 1./D{level+1};
    infNorm = norm(A{level+1}, 'inf');
    
    %---------------
    % update
    %---------------
    if N{level+1}/N{level} >= 0.75
        amgParam.strong_connection = amgParam.strong_connection/4;
    elseif N{level+1}/N{level} >= 0.5 && N{level+1}/N{level} < 0.75
        amgParam.strong_connection = amgParam.strong_connection/2;
    else
        amgParam.strong_connection = amgParam.strong_connection*1.5;
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
    'IL',IL,'IU',IU,'infNorm', infNorm, 'N',N,'max_level',level,...
    'schwarzData',schwarzData,'number_eigen',amgParam.number_eigen);

amgData = amgData(1:level);

% print information
if print_level > 1
   
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
    fprintf('      SA AMG setup costs %f seconds\n', setup_duration);
    fprintf('----------------------------------------------------\n');

end

end

