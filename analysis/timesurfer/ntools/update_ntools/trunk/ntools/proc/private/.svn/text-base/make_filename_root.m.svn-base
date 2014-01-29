function new_root = make_filename_root(old_file)
    [d f e] = fileparts(old_file);
    if(~isempty(e))
        new_root = strip_double_suffix(old_file);
    else
        new_root = old_file;
    end
end