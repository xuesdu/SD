function [ amgData ] = AMG_Setup_Classical( Af, amgParam )
% Setup phase for Classical AMG (CAMG) method
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

%agg_type = amgParam.agg_type;
%agg_radius = amgParam.agg_radius;

MIS_radius = 1;

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
    fprintf('              Calling Classical AMG setup    \n');
    fprintf('----------------------------------------------------\n');
end

setup_start = tic;

while ((level < max_level) && (N{level} > coarsest_size) )
    
    %----------------
    % strong connection
    %----------------
    As = construct_strong_connection( A{level}, amgParam );
    
    % check As
    if (nnz(As) == size(As,1)) 
        break;
    end
    
    %----------------
    % form CF-splitting
    %----------------
    [isC] = coarsening_classical_MIS(As, MIS_radius);
    
    %----------------
    % generate the prolongation
    %----------------
    [ P{level} ] = generate_classical_P_Jacobi(A{level}, isC, MIS_radius);
        
    %----------------
    % generate the restriction
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
    if N{level+1}/N{level} >= 0.5
        amgParam.strong_connection = amgParam.strong_connection*1.2;
    elseif N{level+1}/N{level} >= 0.25 && N{level+1}/N{level} < 0.5
        amgParam.strong_connection = amgParam.strong_connection*1.05;
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
if print_level > 0

    fprintf('----------------------------------------------------\n');
    fprintf('      Classical AMG setup costs %f seconds\n', setup_duration);
    fprintf('----------------------------------------------------\n');

end

end

