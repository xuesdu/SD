function export_MatrixMarket(A, b, mat_file_name, rhs_file_name)

[Nrow,Ncol] = size(A);
Nnz = nnz(A);

[i,j,val] = find(A);
save_data = [i,j,val];      
[~,Index] = sort(save_data(:,1)); 
save_data = save_data(Index,:);

fid = fopen(mat_file_name,'w');
fprintf(fid, '%%MatrixMarket matrix coordinate real general\n');
fprintf(fid, '%d %d %d\n',Nrow,Ncol,Nnz);    
fprintf(fid, '%d %d %0.15e\n', save_data');
fclose(fid);

fid = fopen(rhs_file_name,'w');
fprintf(fid, '%%MatrixMarket matrix array real general\n');
fprintf(fid, '%d %d\n',Nrow, 1);  
fprintf(fid,'%0.15e\n',b');    
fclose(fid);



