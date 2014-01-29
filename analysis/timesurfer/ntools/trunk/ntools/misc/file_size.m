function fsz = file_size(filename)
    fid = fopen(filename, 'r');
    fseek(fid, 0, 'eof');
    fsz = ftell(fid);
    fclose(fid);
end