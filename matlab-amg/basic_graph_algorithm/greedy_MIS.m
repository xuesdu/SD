function [ isM ] = greedy_MIS( A )
% Find maximal indepent set of a given Graph by greedy approach (represented by matrix A)
%
% @ Xiaozhe Hu, Tufts University

% local variables
N = size(A,1);
isM = false(N,1);       % independent set
isS = false(N,1);        % S: selected set

% get permutation
p = amd(A);

% main loop
for i = 1:N
    
    if isS(p(i)) == false
        
        % put this vertex in MIS
        isM(p(i)) = true;
        
        % mark it 
        isS(p(i)) = true;
        
        % mark its neighbors
        [jj,~] = find(A(:,p(i)));
        isS(jj) = true;
        
    end
    
end


end

