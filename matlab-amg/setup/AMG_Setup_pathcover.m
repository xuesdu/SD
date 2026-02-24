function [ amgData ] = AMG_Setup_pathcover( Af, smooth_error, amgParam )
% Setup phase for AMG method based on path covers
%
% @ Xiaozhe Hu, Junyuan Lin (Tufts University) and Ludmil Zikatanov (Penn State)

%----------------
% local variable
%----------------
print_level = amgParam.print_level;

max_level = amgParam.max_level;
coarsest_size = amgParam.coarsest_size;

ILU_level = amgParam.ILU_level;
if ILU_level > 0 
    ILU_setup.type = amgParam.ILU_type;  % nofill, crout, ilutp
    ILU_setup.droptol = amgParam.droptol;
    ILU_setup.milu = amgParam.milu;  % row, col, off
    ILU_setup.udiag = amgParam.udiag;
    ILU_setup.thresh = amgParam.thresh;
end


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
IL = cell(max_level, 1);
IU = cell(max_level, 1);
infNorm = cell(max_level, 1);

%----------------
% finest level
%----------------
A{1} = Af;
N{1} = size(Af,1);
DL{1} = tril(Af);
DU{1} = triu(Af);
D{1} = diag(Af);
Dinv{1} = 1./D{1};
infNorm{1} = norm(Af, 'inf');

if print_level > 0
    fprintf('----------------------------------------------------\n');
    fprintf('     Calling AMG setup based on path covering    \n');
    fprintf('----------------------------------------------------\n');
end
    
setup_start = tic;

%----------------
% get paths based on smooth_error
%----------------
% compute affinity
%[ A_affinity ] = compute_affinity(smooth_error, Af, 0.01);
[ A_affinity ] = compute_similarity( smooth_error, Af, 2 );
%[ A_affinity ] = compute_weights(smooth_error, Af, 8 );
%[ A_affinity ] = compute_similarity_knn(smooth_error, 1 );
As = -A_affinity;

% compress As
[row, col, val] = find(triu(As,1));
[~,idx] = sort(val);
start = floor(length(val)*0.0)+1;
As = sparse(row(idx(start:end)), col(idx(start:end)), val(idx(start:end)), N{1}, N{1});
As = As+As';

% generate the cover
[ cover ] = genCover1( As );


%----------------
% main loop
%----------------
while ((level < max_level) && (N{level} > max(coarsest_size, 0)) )
    
    %----------------
    % form aggregation along the paths
    %----------------
    [ aggregation,num_agg, iso_edges ] = pair_match( cover, A_affinity );
    [ cover ] = coarsen_pathcover(cover, aggregation);
    
    %----------------
    % generate prolongation
    %----------------
    if  (level == 1)
        [ P{level} ] = generate_unsmoothed_P_smooth_error( aggregation, num_agg, smooth_error);
    else
        [ P{level} ] = generate_unsmoothed_P(aggregation, num_agg);
    end
        
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
    % update affinity
    %----------------
    A_affinity = R{level}*A_affinity*P{level};
    
    %----------------
    % extra information
    %----------------
    if level <= ILU_level
        [IL{level}, IU{level}] = ilu(A{level}, ILU_setup); 
    end
        
    DL{level+1} = tril(A{level+1});
    DU{level+1} = triu(A{level+1});
    D{level+1} = diag(A{level+1});
    Dinv{level+1} = 1./D{level+1};
    infNorm = norm(A{level+1}, 'inf');
    
    %---------------
    % update
    %---------------
    if N{level+1}/N{level} >= 0.9
        amgParam.strong_connection = amgParam.strong_connection/8;
    elseif N{level}/N{level} >= 0.5 && N{level+1}/N{level} < 0.9
        amgParam.strong_connectionn = amgParam.strong_connection/4;
    end
    level = level+1;
        
end

setup_duration = toc(setup_start);

% construct the data structure
amgData = struct('P',P,'R',R,'A',A,'DL',DL,'DU',DU,'D',D,'Dinv',Dinv,...
    'IL',IL,'IU',IU,'infNorm', infNorm, 'N',N,'max_level',level);

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

if print_level > 0
    fprintf('----------------------------------------------------\n');
    fprintf('    path cover AMG setup costs %f seconds\n', setup_duration);
    fprintf('----------------------------------------------------\n');
end

end

