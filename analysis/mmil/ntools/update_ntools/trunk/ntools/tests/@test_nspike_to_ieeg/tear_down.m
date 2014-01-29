function self = tear_down(self)
    fclose(self.fid_nspike);
    %fclose(fid_nspike);
    self.orig_data = [];
end
