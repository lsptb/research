% predefine loadgrad, loadeeg, writelog, outpath
subject   = 's1';
loadflag  = loadgrad || loadeeg;
if writelog
  logfile = sprintf('%s/%s.log',outpath,date);
  fid     = fopen(logfile,'a');
  fprintf(fid,'---------------------------------\n%s\nSlow oscillation analysis\n---------------------------------\n',datestr(now));
else
  fid = 1;
end
if loadflag
  findex    = 1:5;
    if strcmp(subject,'s1')
      badlabels = {'C1','Cz','C6','CPz','CP4'};
      fiffiles  = {...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_2_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_3_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_4_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_5_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_6_raw.fif'};% ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_DC_s1_7_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_DC_s1_8_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_9_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_10_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_11_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_12_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_13_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_14_raw.fif'};
    elseif strcmp(subject,'s2')
      badlabels = {'MEG 2112','EEG 018','EEG 029','EEG 033','EEG 037','EEG 042'};
      fiffiles  = {...        
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_1_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_2_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_3_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_4_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_5_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_6_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_7_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_8_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_9_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_DC_s2_12_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_DC_s2_13_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_10_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_11_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_14_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_15_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_16_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_17_raw.fif' ...      
    elseif strcmp(subject,'s8')
      fiffiles  = {...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_1_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_2_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_3_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_4_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_5_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_6_nb01_060808.fif' ...
    %     '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/emptyroom_nb01_060808.fif' ...
        };  
    else
      error('Files have not been specified for %s',subject);
    end      
  % read fif files, convert data to timesurfer format, and save mat files
  matfiles = {};
  chantype = {'grad1','grad2'};
  for f = 1:length(fiffiles)
    fif = fiffiles{f};
    [fpath,fname,fext]  = fileparts(fif);
    outfile             = sprintf('%s/matfiles/%s_grad.mat',outpath,fname);
    matfiles{end+1}     = outfile;
    if exist(outfile,'file') % never overwrite (param independent)
      fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
      continue
    else
      fprintf(fid,'Reading FIF file: %s\n',fif);
    end
    data = ts_loadfif(fif,chantype,'epochs');
    fprintf(fid,'Saving MAT file: %s\n',outfile);
    save(outfile,'data');
    clear fif data
  end
  if loadgrad
    clear data
    % read mat files and combine data
    fprintf(fid,'Loading MAT files:\n');
    for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
    data   = SO_combine_matfiles(matfiles(findex));
  end
  if loadeeg
    clear eeg
    for f = 1:length(fiffiles)
      fif = fiffiles{f};
      [fpath,fname,fext]  = fileparts(fif);
      outfile             = sprintf('%s/matfiles/%s_eeg.mat',outpath,fname);
      matfiles{end+1}     = outfile;
      if exist(outfile,'file') % never overwrite (param independent)
        fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
        continue
      else
        fprintf(fid,'Reading FIF file: %s\n',fif);
      end
      data = ts_loadfif(fif,{'eeg'},'epochs');
      fprintf(fid,'Saving MAT file: %s\n',outfile);
      save(outfile,'data');
      clear fif data
    end    
    findex = length(fiffiles) + findex;
    for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
    eeg    = SO_combine_matfiles(matfiles(findex));    
  end
end
% peaks for grad and/or eeg
load(eventfile); % events, peaks
% example: s1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat

clear f fext fif findex fname fpath k outfile