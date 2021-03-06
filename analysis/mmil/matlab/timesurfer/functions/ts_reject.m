function [data reject_data] = ts_reject(data,varargin)
% data/outdata: avg_data, epoch_data, timefreq_data
% note: needs to work for both parameter names from iEEG and fif streams
%
% Flags [0|1]:
% threshold_flag
% visualreject_flag
% ICA_auto_flag
% ICA_manual_flag
%
%  'reject_mag'  - automatic rejection threshold for magnetometer channels (fT)
%     if 0, rejection based on magnetometers is disabled
%     { default: 6000 }
%  'reject_grad' - automatic rejection threshold for gradiometer channels (fT/cm)
%     if 0, rejection based on gradiometers is disabled
%     { default: 3000 }
%  'reject_eeg' - automatic rejection threshold for eeg channels (uV)
%     if 0, rejection based on eeg is disabled
%     { default: 0 }
%  'reject_eog' - automatic rejection threshold for eog channel (uV)
%     if 0, rejection based on eog is disabled
% 
% Default: return reject_data (w/o saving) and data (badtrials marked bad
% and minus the removed ICs)

% data  = ts_checkdata_header(data);
parms = mmil_args2parms(varargin,...
						{'save_reject_flag'   ,0,{0,1},...
             'keepbadtrials_flag' ,1,{0,1},...
             'visualreject_flag'  ,0,{0,1},...
             'threshold_flag'     ,0,{0,1},...
             'ICA_auto_flag'      ,0,{0,1},...
             'ICA_manual_flag'    ,0,{0,1},...
             'ICA_ref_chan'       ,'EOG061',[],...
             'ICA_chantype'       ,'all',[],...
             'ICA_maxsteps'       ,20,[],...
             'ICA_ntrial'         ,5,[],...
             'ICA_ncomponents'    ,80,[],...
             'ICA_plottype'       ,'activations',{'activations','alltrials'},...
             'ICA_rescale_flag'   ,1,{0,1},...
             'ICA_sorttrials'     ,0,{0,1},...       
             'reject_method'      ,'summary',[],...
             'reject_metric'      ,'var',[],...     
             'rejectfile'         ,[],[],...
             'reject_data'        ,[],[],...
             'reject_mag'         ,6000,[],...
             'reject_grad'        ,3000,[],...
             'reject_eeg'         ,0,[],...
             'reject_eog'         ,200,[],...      
             'reject_ieeg'        ,[],[],...
             'prescale_mag'       ,10^15,[],...  % to fT
             'prescale_grad'      ,10^13,[],... % to fT/cm
             'prescale_eeg'       ,10^6,[],...   % to uV
             'prescale_eog'       ,10^6,[],...   % to uV
             'prefix','proc'      ,[],...
             'events'             ,[],[],...
						 'conditions'         ,[],[],...
             'badchanfile'        ,[],[],...
             'badchans'           ,[],[],...
             'toilim'             ,[],[],...
             'verbose'            ,1,{0,1},...
             'logfile'            ,[],[],...
             'logfid'             ,[1],[], ...  
             'filename'           ,[],[],...
             'write_reject_log'   ,0,{0,1},...
						},false);

parms = backcompatible(parms,varargin{:});      
data  = ts_checkdata_header(data,'events',parms.events);
if exist(parms.badchanfile,'file') || ~isempty(parms.badchans) || ~isempty(parms.toilim)
  data = ts_data_selection(data,'rejectfile',parms.rejectfile,'reject_data',parms.reject_data,'badchanfile',parms.badchanfile,'badchans',parms.badchans,'toilim',parms.toilim);
end
[datatype,datafield,dataparam] = ts_object_info(data);

nchan  = data.num_sensors;
ncond  = length(data.(datafield));
[allbadtrials{1:ncond}] = deal([]);
[allbadsamps{1:ncond}]  = deal([]);
allbadchans             = [];
%% Automatic Threshold Rejection
if parms.threshold_flag
  % determine which channels are which type
  typestring = {data.sensor_info.typestring};
  mag_i   = strcmp ('mag' ,typestring);
  grad_i  = strncmp('grad',typestring,4);
  eeg_i   = strcmp ('eeg' ,typestring);
  eog_i   = strcmp ('eog' ,typestring);
  % generate rejection thresholds for each channel
  reject_thresh         = zeros(nchan,1);
  reject_thresh(mag_i)  = parms.reject_mag/parms.prescale_mag;
  reject_thresh(grad_i) = parms.reject_grad/parms.prescale_grad;
  reject_thresh(eeg_i)  = parms.reject_eeg/parms.prescale_eeg;
  reject_thresh(eog_i)  = parms.reject_eog/parms.prescale_eog;
  reject_thresh(reject_thresh<=0) = inf;
  for c = 1:ncond
    ntrl      = data.(datafield)(c).num_trials;
    nsamp     = length(data.(datafield)(c).time);
    tmpthresh = reject_thresh * ones(1,nsamp);
    rej_idx   = cellfun(@(x)find(abs(data.(datafield)(c).data(:,:,x))>tmpthresh),num2cell(1:ntrl),'uniformoutput',false);
    badtrl_ix = find(~cellfun(@isempty,rej_idx));
    badchan_ix= [];
    for b = 1:length(badtrl_ix)
      trl = badtrl_ix(b);
      [chan_ix,samp_ix] = ind2sub([nchan nsamp],rej_idx{trl});
      badchan_ix(b)     = chan_ix(1);
    end
    allbadtrials{c} = [allbadtrials{c} badtrl_ix];
    allbadsamps{c}  = [allbadsamps{c}  data.(datafield)(c).trial_info.latency(badtrl_ix)];
    badchantypes    = typestring(badchan_ix);
    tmp       = data.(datafield)(c).num_rejects;
    tmp.mag   = tmp.mag  + length(find(strmatch('mag' ,badchantypes)));
    tmp.grad  = tmp.grad + length(find(strmatch('grad',badchantypes)));
    tmp.eeg   = tmp.eeg  + length(find(strmatch('eeg' ,badchantypes)));
    tmp.eog   = tmp.eog  + length(find(strmatch('eog' ,badchantypes)));
    data.(datafield)(c).num_rejects = tmp;
    clear tmp
    data.(datafield)(c).trial_info.badtrial(badtrl_ix) = 1;
    % TODO: add reject log (take from Sanja's ts_process_fif_data)
  end
end

%% Manual Visual Rejection

% ts_browseraw()
% ts_browseraw_graph()

% ts_visual_reject
if parms.visualreject_flag
  remove_bad_trials; % remove bad trials before visual rejection
  ix = find(cellfun(@(x)isequal(x,'method'),varargin)); if any(ix), varargin([ix ix+1])=[]; end
  ix = find(cellfun(@(x)isequal(x,'metric'),varargin)); if any(ix), varargin([ix ix+1])=[]; end
  [data,badchans,badtrials] = ts_visual_reject(data,varargin{:},'method',parms.reject_method,'metric',parms.reject_metric);
  allbadchans  = [allbadchans badchans];
  for k = 1:ncond
    if isempty(badtrials{k}), continue; end
    allbadtrials{k} = [allbadtrials{k} badtrials{k}];
    allbadsamps{k}  = [allbadsamps{k}  data.(datafield)(k).trial_info.latency(badtrials{k})];
    data.(datafield)(k).trial_info.number(badtrials{k})        = [];
    data.(datafield)(k).trial_info.latency(badtrials{k})       = [];
    data.(datafield)(k).trial_info.badtrial(badtrials{k})      = [];
    data.(datafield)(k).trial_info.event_code(badtrials{k})    = [];
    data.(datafield)(k).trial_info.duration(badtrials{k})      = [];
    data.(datafield)(k).trial_info.datafile(badtrials{k})      = [];
    data.(datafield)(k).trial_info.events_fnames(badtrials{k}) = [];    
  end
end

%% Automtic ICA
if parms.ICA_auto_flag
  % perform automatic ICA rejection for blinks
  data = ts_autoICA(data,'ICA_ref_chan',parms.ICA_ref_chan,'chantype',...
                    parms.ICA_chantype,'rescale',parms.ICA_rescale_flag);
end

%% Manual ICA
if parms.ICA_manual_flag
  % perform manual ICA rejection for blinks and/or EKG
  data = ts_manualICA(data,'maxsteps',parms.ICA_maxsteps,'ntrial',parms.ICA_ntrial,...
                'ncomponents',parms.ICA_ncomponents,'chantype',parms.ICA_chantype,'plottype',parms.ICA_plottype,...
                'rescale',parms.ICA_rescale_flag,'sorttrials',parms.ICA_sorttrials,'allconditions',1);
end

%% Create and save Reject Data (*.mat, *.log)
if 1
  reject_data = [];
  reject_data.badchans      = allbadchans;
  reject_data.badchanlabels = {data.sensor_info(allbadchans).label};
  reject_data.badtrials= [];
  for i = 1:length(data.(datafield))
    reject_data.badtrials{i}  = allbadtrials{i};
    reject_data.badsamples{i} = allbadsamps{i}; %data.(datafield)(i).trial_info.latency(allbadsamps{i});
    reject_data.event_code(i) = data.(datafield)(i).event_code;
  end
else
  reject_data = [];
end

if parms.save_reject_flag
  % Save badchannel and trials as matlab file
  tmp        = rmfield(data,'epochs');
  tmp.epochs = rmfield(data.epochs,'data');
  save_reject(reject_data,parms,tmp);  
  clear tmp
end

if ~parms.keepbadtrials_flag
  remove_bad_trials;
end

% if parms.write_reject_log
%   write_reject_log;
% end

  function remove_bad_trials
    for c = 1:ncond
      % remove rejects from data
      idx = find(data.(datafield)(c).trial_info.badtrial == 1);
      if isempty(idx) || all(idx==0), continue; end
      data.(datafield)(c).trial_info.number(idx)        = [];
      data.(datafield)(c).trial_info.latency(idx)       = [];
      data.(datafield)(c).trial_info.badtrial(idx)      = [];
      data.(datafield)(c).trial_info.event_code(idx)    = [];
      data.(datafield)(c).trial_info.duration(idx)      = [];
      data.(datafield)(c).trial_info.datafile(idx)      = [];
      data.(datafield)(c).trial_info.events_fnames(idx) = [];
      data.(datafield)(c).data(:,:,idx) = [];
      data.(datafield)(c).num_trials = data.(datafield)(c).num_trials - length(idx);
    end
    rmix                    = find([data.(datafield).num_trials]==0);
    ncond                   = ncond - length(rmix);
    allbadtrials(rmix)      = [];
    allbadsamps(rmix)       = [];
    data.(datafield)(rmix)  = [];
  end

%   function write_reject_log
%     % where to save?
%     if ~isempty(parms.rootoutdir)
%       outdir = parms.rootoutdir;
%     elseif ~isempty(parms.evntfile)
%       outdir = fileparts(parms.evntfile);
%     else
%       outdir = pwd;
%     end
%     % what to name it?
%     [jnk, name, ext] = fileparts(parms.datafile);
%     if isempty(ext) , ext  = 'TimeSurfer'; end
%     if isempty(name), name = sprintf('%s',parms.prefix); end
%     logfile = sprintf('%s/%s_events.log',outdir,name);
%     % write the log
%     fid     = fopen(logfile,'wt');
% %     [~,b]   = unix('echo $HOST');
% %     fprintf(fid,'--------------------------------------------------------------------\n');
% %     fprintf(fid,'Time: %s\n',datestr(now));
% %     fprintf(fid,'Host: %s',b);
% %     fprintf(fid,'Epoching %s data:\n%s\n',ext,parms.datafile);
% %     fprintf(fid,'Number channels: %i \n',data.num_sensors);
% %     fprintf(fid,'Sampling frequency: %5.6g Hz\n',Fs);
% %     fprintf(fid,'--------------------------------------------------------------------\n');
%     fprintf(fid, '%-6s\t %-10s\t %-10s\t %-12s\t %-6s\t %-20s\n','Trial','Event','Latency','Time (sec)','Type','Comment');
%     for k = 1:nevnt
%         fprintf(fid,'%-6g\t %-10g\t %-10g\t %-12.8g\t %-6s\t %-20s\n',k,cond(k),samp(k),samp(k)/Fs,'ok','no rejection');
%     end
% %     fprintf(fid,'--------------------------------------------------------------------\n');
%     fclose(fid);
%   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_reject(reject_data,parms,hdr)
  % save reject_data
  % Save reject data
  if ~isempty(reject_data)
    % save info to mat file
    if isempty(parms.filename)
      args           = mmil_parms2args(parms); 
      parms.filename = ts_create_filename('ts_reject',args{:});
    end
    if iscell(parms.filename)
      filename = parms.filename{1};
    else
      filename = parms.filename;
    end
    filename = strrep(filename,'epoch_data_rej','reject_data');
    mmil_logstr(parms,'%s: saving reject data: %s\n',mfilename,filename);
    save(filename,'reject_data');

    % save info to text file (copied from Rajan's code in ts_iEEG_ProcessEEG)
    filename = strrep(filename,'.mat','.txt');
    fid = fopen(filename,'w+');
    badchans = find([hdr.sensor_info.badchan]);  
    fprintf(fid,'Rejected Channels\n\n');
    if ~isempty(badchans)
      for j = 1:length(badchans)
        fprintf(fid,'%s\n',hdr.sensor_info(badchans(j)).label);
      end
    else
      fprintf(fid,'No bad channels.\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Trial Info\n\n');
    fprintf(fid,'         \t Good\t    Rejected     \t Bad\n');
    fprintf(fid,'Condition\tTrials\tEEG File\tManual\tTrials\n\n');
    for j = 1:length(hdr.epochs)
      if j > length(reject_data.badtrials), break; end
     fprintf(fid,'%-9s\t%-6s\t%-8s\t%-6s\t%s\n',num2str(j),num2str(hdr.epochs(j).num_trials),...
                                                       num2str(hdr.epochs(j).num_rejects.eeg),...
                                                       num2str(hdr.epochs(j).num_rejects.manual),...
                                                       num2str(reject_data.badtrials{j}));
    end
    fclose(fid);  
  end
end
%%
function parms = backcompatible(parms,varargin)

opt = mmil_args2parms(varargin,{...
     'ICA_ntrials','',[],...   
     'method',[],[],...
     'metric',[],[],...
     'visualreject',[],[],...
     },false);

if ~isempty(opt.ICA_ntrials)
  parms.ICA_ntrial = opt.ICA_ntrials;
end
if ~isempty(opt.method)
  parms.reject_method = opt.method;
end
if ~isempty(opt.metric)
  parms.reject_metric = opt.metric;
end
if isempty(parms.visualreject_flag)
  parms.visualreject_flag = opt.visualreject;
end
end