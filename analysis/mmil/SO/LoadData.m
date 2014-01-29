% Purpose: this script loads the raw data from mat files with data in
% TimeSurfer format (continuous data in an epoch_data structure with
% one trial), concatenates, preprocesses and runs ICA to remove EKG.  The script
% will start with the most recently saved files (ex. the post-ICA
% epoch_data will be loaded directly if the MAT file already exists).
% Note: use ConvertFIF2MAT.m before calling this function.
%
% Created by Jason Sherfey on 27-Jul-2010
% Multimodal Imaging Laboratory
% Department of Radiology, UCSD

outfile = sprintf('matfiles/proc_%s_epoch_data_ICA.mat',params.chantype);
if exist(outfile,'file')
  fprintf(fid,'Loading post-ICA continuous %s sleep data: %s\n',params.chantype,outfile);
  load(outfile,'epoch_data');
  data = epoch_data;
  clear epoch_data  
  if ~isempty(params.toilim)
    fprintf(fid,'Selecting %g-%gsec from %g-%gsec\n',params.toilim,data.epochs.time(1),data.epochs.time(end));
    data = ts_data_selection(data,'toilim',params.toilim);
  end  
else
  file = sprintf('matfiles/proc_%s_epoch_data_1.mat',params.chantype);
  if exist(file,'file')
    fprintf(fid,'Loading pre-ICA continuous sleep data: %s\n',file);
    load(file,'epoch_data'); clear file
  else
    % read mat files and concatenate data
    fprintf(fid,'Loading MAT files:\n');
    for k  = 1:length(params.matfile_index)
      fprintf(fid,'%s\n',params.datafiles{params.matfile_index(k)});
    end
    if length(params.matfile_index) > 10 && strcmp(params.chantype,'grad')
      epoch_data   = SO_combine_matfiles(params.datafiles(params.matfile_index),params,fid,1,params.chantype);
      if ~isempty(params.toilim)
        % each matfile will have been preprocessed indiv by SO_combine_matfiles
        fprintf(fid,'Selecting times %g-%gsec.\n',params.toilim);
        epoch_data = ts_data_selection(epoch_data,'toilim',params.toilim);
      end
      fprintf(fid,'Sampling rate = %gHz.\n',epoch_data.sfreq);
    else
      data   = SO_combine_matfiles(params.datafiles(params.matfile_index),[],1,1,params.chantype);
      if ~isempty(params.toilim)
        fprintf(fid,'Selecting times %g-%gsec.\n',params.toilim);
        data = ts_data_selection(data,'toilim',params.toilim);
      end
      fprintf(fid,'Sampling rate = %gHz.\n',data.sfreq);
      fprintf(fid,'Concatenating and preprocessing data before ICA.\n');
      epoch_data = ts_preproc(data,   'dsfact',     params.ICA_dsfact,      ...
                                      'bpfilter',   'yes',  'bpfreq',   params.ICA_bpfreq, 'bandpass_detrend_flag',0,...
                                      'notch_flag', params.ICA_notchflag,      ...
                                      'blc',        params.ICA_blcflag,  'blcwindow',params.ICA_blcwindow  );        
    end
    fprintf(fid,'\tConcatenated time interval: %g-%gsec.\n',[epoch_data.epochs.time(1) epoch_data.epochs.time(end)]);
    fprintf(fid,'\tDownsample factor = %g.\n',params.ICA_dsfact);
    fprintf(fid,'\tEffective sampling rate = %gHz.\n',epoch_data.sfreq);
    fprintf(fid,'\tNot removing linear fit before filtering.\n');
    fprintf(fid,'\t%g-%gHz zero phase shift frequency domain bandpass filter.\n',params.ICA_bpfreq);
    if params.ICA_notchflag,             fprintf(fid,'\tApply 60Hz notch filter.\n'); end
    if strcmp(params.ICA_blcflag,'yes'), fprintf(fid,'\tRemove mean: [%g %g] sec.\n',params.ICA_blcwindow); end
    clear matfiles fif fpath fname fext
    % Prepare data for ICA removal of EKG artifact
    % eliminate (outer) edge effects
    epoch_data = ts_data_selection(epoch_data,'toilim',[round(epoch_data.epochs.time(1)+5) round(epoch_data.epochs.time(end)-5)]);
    matfile    = sprintf('%s/matfiles/proc_%s_epoch_data_1.mat',params.SubjDir,params.chantype); 
    if strcmp(params.chantype,'eeg') && isfield(params,'eeg_index') && ~isempty(params.eeg_index)
      fprintf(fid,'Rearranging EEG channels before saving data structure.\n\tShuffle: [%s]\n',num2str(params.eeg_index));
      epoch_data.sensor_info = epoch_data.sensor_info(params.eeg_index);
      epoch_data.num_sensors = length(epoch_data.sensor_info);
      epoch_data.epochs.data = epoch_data.epochs.data(params.eeg_index,:);
    end
    fprintf(fid,'Saving matfile before running ICA: %s\n',matfile);
    save(matfile,'epoch_data','-v7.3');
  end
  % ICA removal of EKG artifact
  fprintf(fid,'Running ICA for times: %g-%gsec.\n',epoch_data.epochs.time(1),epoch_data.epochs.time(end));
  data       = ts_manualICA(epoch_data,'maxsteps',params.ICA_maxsteps,'logfid',fid,...
    'ntrial',params.ICA_ntrial/epoch_data.epochs.time(end),'prefix',sprintf('proc_%s',params.chantype));
  fprintf(fid,'New ICA epoch_data file: %s/matfiles/proc_%s_epoch_data_ICA.mat\n',params.SubjDir,params.chantype);
  clear epoch_data matfile
end
clear outfile