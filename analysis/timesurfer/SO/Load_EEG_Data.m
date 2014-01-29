addpath(genpath('/usr/pubsw/packages/mne-2.6.1_64/mne'));
params = SO_params(1);
params.monotonic = 0;
fid = 1;
typestr   = 'grad2'; % 'eeg';
nogap     = 0;
begsamp   = 0;
matfiles  = {};
fprintf(fid,'Processing %s data\n',typestr);%[params.chantype{:}]);
for f = 1:length(params.datafiles)
  if ~(f>=min(params.matfile_index) && f<=max(params.matfile_index)), continue; end %ismember(f,params.matfile_index), continue; end
  fif = params.datafiles{f};
  [fpath,fname,fext]  = fileparts(fif);
  outfile             = sprintf('%s/matfiles/%s_%s.mat',params.SubjDir,fname,typestr);
  matfiles{end+1}     = outfile;
  if exist(outfile,'file') % never overwrite (param independent)
%       fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
    continue
  else
    fprintf(fid,'Reading FIF file: %s\n',fif);
  end
  % ...
%     data = ts_loadfif(fif,params.chantype,'epochs');

  endsamp = [];
  hdr     = ts_read_fif_header(fif,0);
  for k   = 1:length(hdr.sensors.label)
    sens(k).label         = hdr.sensors.label{k};
    sens(k).typestring    = hdr.sensors.typestring{k};
    sens(k).type          = hdr.sensors.type(k);
    sens(k).kind          = hdr.sensors.kind(k);
    sens(k).badchan       = 0;
    sens(k).lognum        = hdr.sensors.lognum(k);
    sens(k).loc           = hdr.sensors.loc{k};
  end
  chans = strmatch(typestr,{sens.typestring});
  sens  = sens(chans);
  hdr   = fiff_setup_read_raw(fif);
  if isempty(endsamp)
    % set endsamp to total number of samples
    endsamp = hdr.last_samp;
  end
  [x,t] = fiff_read_raw_segment(hdr,begsamp,endsamp,chans);
  data  = ts_matrix2epoch(x,'continuous',1,'sens',sens,'time',t);
  clear hdr x t sens chans
  fprintf(fid,'Saving MAT file: %s\n',outfile);
  save(outfile,'data');
  clear fif data
end
% read mat files and concatenate data
fprintf(fid,'Loading MAT files:\n');
for k  = 1:length(params.matfile_index)
  fprintf(fid,'%s\n',matfiles{params.matfile_index(k)});
end
data   = SO_combine_matfiles(matfiles(params.matfile_index),[],fid,nogap);
if ~isempty(params.toilim)
  data = ts_data_selection(data,'toilim',params.toilim);
end

fprintf(fid,'Preprocessing data before ICA\n');
if f > 10
  divind     = round(data.num_sensors/2);
  tmpdat     = ts_preproc(ts_data_selection(data,'channels',1:divind),'dsfact',     params.ICA_dsfact,      ...
                                'precision',  'single',...
                                'bpfilter',   'yes',  'bpfreq',   params.ICA_bpfreq, 'bandpass_detrend_flag',0,...
                                'notch_flag', params.ICA_notchflag,      ...
                                'blc',        params.ICA_blcflag,  'blcwindow',params.ICA_blcwindow  );        
  matdat(1:divind,:) = tmpdat.epochs.data;
  clear tmpdat
  tmpdat     = ts_preproc(ts_data_selection(data,'channels',divind+1:data.num_sensors),'dsfact',     params.ICA_dsfact,      ...
                                'precision',  'single',...
                                'bpfilter',   'yes',  'bpfreq',   params.ICA_bpfreq, 'bandpass_detrend_flag',0,...
                                'notch_flag', params.ICA_notchflag,      ...
                                'blc',        params.ICA_blcflag,  'blcwindow',params.ICA_blcwindow  );        
  data.epochs = rmfield(data.epochs,'data');
  epoch_data  = data;
  epoch_data.sfreq        = tmpdat.sfreq;
  epoch_data.epochs.time  = tmpdat.epochs.time;
  matdat(divind+1:epoch_data.num_sensors,:) = tmpdat.epochs.data;
  clear tmpdat               
  epoch_data.epochs.data = matdat;
  clear matdat
else  
  epoch_data = ts_preproc(data,   'dsfact',     params.ICA_dsfact,      ...
                                'precision',  'single',...
                                'bpfilter',   'yes',  'bpfreq',   params.ICA_bpfreq, 'bandpass_detrend_flag',0,...
                                'notch_flag', params.ICA_notchflag,      ...
                                'blc',        params.ICA_blcflag,  'blcwindow',params.ICA_blcwindow  );        
end
clear matfiles fif fpath fname fext
% Prepare data for ICA removal of EKG artifact
% eliminate (outer) edge effects
epoch_data = ts_data_selection(epoch_data,'toilim',[round(epoch_data.epochs.time(1)+5) round(epoch_data.epochs.time(end)-5)]);
matfile    = sprintf('%s/matfiles/proc_%s_epoch_data_1.mat',params.SubjDir,typestr); 
fprintf(fid,'Saving matfile: %s\n',matfile);
save(matfile,'epoch_data','-v7.3');

data = epoch_data; clear epoch_data
params.monotonic = 0;
args              = mmil_parms2args(params);
detections        = SO_detection(data,args{:});

t                 = data.epochs.time;
events  = [];
  % The events structure can be used by visualizer to mark detections
for k   = 1:length(detections)
  events(k).label = data.sensor_info(k).label;
  events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
  events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
end

save(sprintf('SO_detections_%s_with_gaps.mat',typestr),'events','detections','params');

%% remove gaps
% load proc_eeg_epoch_data_1 % epoch_data
load(sprintf('proc_%s_epoch_data_1.mat',typestr));
epoch_data.epochs.time = ([0:length(epoch_data.epochs.time)-1]/epoch_data.sfreq)+5;
tau = epoch_data.epochs.duration;
mat = epoch_data.epochs.matfiles;
for k = 1:length(mat)
  tmp = load(mat{k});
  beg(k) = tmp.data.epochs.time(1);
  clear tmp
end
% beg(1)  = 5; % 5 sec of padding removed above
tau     = tau - beg;
for k   = 1:length(events)
  T     = events(k).time;
  for j = 1:length(beg)
    ind = find(events(k).time>=beg(end-j+1));
    T(ind) = T(ind) - beg(end-j+1);
  end
  events(k).time = T;
end
save(sprintf('SO_detections_%s_without_gaps.mat',typestr),'events','detections','params');

%% temporary: copied from combine_meg_eeg
%% combine grad and eeg events
load('SO_detections_eeg_without_gaps.mat','events'); tmpevents = events;
load('SO_detections_grad2_without_gaps.mat','events');
events(end+1:end+length(tmpevents)) = tmpevents;
clear tmpevents
save('SO_events_meg_eeg_without_gaps.mat','events','params'); 

%% combine grad and eeg data for visualization
load proc_grad2_epoch_data_1; data=epoch_data;
load proc_eeg_epoch_data_1; eeg=epoch_data;
clear epoch_data

ng  = data.num_sensors;
ne  = eeg.num_sensors;
tg  = data.epochs.time;
te  = eeg.epochs.time;
  % eeg time vector is wrong (need to investigate)
  % use grad time vector for now
  % note: there's an 85ms discrepancy...

data.sensor_info(ng+1:ng+ne) = eeg.sensor_info;
data.num_sensors = length(data.sensor_info);

ta  = 1;
tb  = min(length(tg),length(te));

data.epochs.time = data.epochs.time(ta:tb);
data.epochs.data = data.epochs.data(:,ta:tb);
data.epochs.data(ng+1:ng+ne,:) = eeg.epochs.data(:,ta:tb);
data.epochs.time = ([0:length(data.epochs.time)-1]/data.sfreq)+5;
visualizer(data); % load combined events structure (SO_events_meg_eeg_without_gaps.mat)

%% 
load('SO_detections_eeg_with_gaps.mat','events'); tmpevents = events;
load('SO_detections_grad2_with_gaps.mat','events');
events(end+1:end+length(tmpevents)) = tmpevents;
clear tmpevents
save('tmpevents.mat','events');
visualizer(data); % load tmpevents.mat
