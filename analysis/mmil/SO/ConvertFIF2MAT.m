scriptpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/Cancellation';
cwd = pwd;
cd(scriptpath);

tic
outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/Cancellation';
Raw_FIF_files;
% fiffiles = {'/home/halgdev/data/MEG_UCSD/SL_nd01_060726/wake_nd01_060726.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_0_nd01_060726.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_1_nd01_060726.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/awake_tb01_060801.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_1_tb01_060801.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/awake_ma01_060729.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_0_ma01_060729.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_1_ma01_060729.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_5_ma01_060729.fif',...
%             '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_6_ma01_060729.fif'};          
for f = 1:length(fiffiles)
  fiffile       = fiffiles{f};
  [jnk,outfile] = fileparts(fiffile); 
  outfile       = sprintf('%s/%s.mat',outpath,outfile);
  if exist(outfile,'file'), continue; end
  data = ts_MNE_loadfif(fiffile);
  save(outfile,'data','-v7.3');
  clear data
end

toc

cd(cwd); 