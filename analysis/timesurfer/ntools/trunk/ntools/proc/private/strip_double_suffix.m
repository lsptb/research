function newfile = strip_double_suffix(oldfile)
    [d f e] = fileparts(oldfile);
    [dummy ff] = fileparts(f);
    newfile = fullfile(d, ff);
end