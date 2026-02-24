function export(A, b)
%function export(A)
%function export(b)

%A = sparse(J);
 [Nrow,Ncol] = size(A);
 Nnz = nnz(A);
% 
 [i,j,val] = find(A);
 save_data = [i-1,j-1,val];   % C friendly index
% %save_data = [i,j,val];      % Nature index
 [temp,Index] = sort(save_data(:,1)); % sort by row, which fits read_IJ_matrx in MSC package
 save_data = save_data(Index,:);
% 
% %fid = fopen([[['./Test/ijmatJ_',num2str(n)],'_',num2str(m)],'_',num2str(Nrow)],'w');
% %fid = fopen('./schur/model2_schur_mat.dat','w'); 
% %fid = fopen('./ml/model3_homo_aij.dat','w');
 fid = fopen('./P.dat','w');
% %------------------------------------------------
% %------------------------------------------------
 fprintf(fid,'%d %d %d\n',Nrow,Ncol,Nnz);    % msc
% %------------------------------------------------
% %fprintf(fid, '%d\n', Nrow); %ml
% %fprintf(fid, '%d\n', Nnz); %ml
% %fprintf(fid, '%d\n', 0); %ml
% %------------------------------------------------
 fprintf(fid, '%d %d %0.15e\n', save_data');  
 fclose(fid);

%fid = fopen([[['./Test/ijrhs_',num2str(n)],'_',num2str(m)],'_',num2str(Nrow)],'w');
%fid = fopen('./schur/model2_schur_rhs.dat','w'); 
%fid = fopen('./ml/model3_homo_rhs.dat','w');
fid = fopen('./d.dat','w');
%-------------------------------------------------
%-------------------------------------------------
fprintf(fid, '%d\n',Nrow);  
fprintf(fid,'%0.15e\n',b');    % msc
fclose(fid);



