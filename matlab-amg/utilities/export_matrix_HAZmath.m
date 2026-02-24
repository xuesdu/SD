function export_matrix_HAZmath(A, filename)

    [Nrow,Ncol] = size(A);
    Nnz = nnz(A);

    [i,j,val] = find(A);
    save_data = [i-1,j-1,val];   % C friendly index

    [~,Index] = sort(save_data(:,1)); % sort by row, which fits HAZmath package
    save_data = save_data(Index,:);
    
    fid = fopen(filename, 'w');
    
    % %------------------------------------------------
    % %------------------------------------------------
    fprintf(fid,'%d %d %d\n',Nrow,Ncol,Nnz);    % size
    % %------------------------------------------------
    fprintf(fid, '%d %d %0.15e\n', save_data');  
    
    fclose(fid);




