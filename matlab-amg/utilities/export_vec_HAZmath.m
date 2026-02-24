function export_vec_HAZmath(b, filename)

    Nrow = length(b);

    fid = fopen(filename,'w');
    %-------------------------------------------------
    %-------------------------------------------------
    fprintf(fid, '%d\n',Nrow);  
    fprintf(fid,'%0.15e\n',b');

fclose(fid);



