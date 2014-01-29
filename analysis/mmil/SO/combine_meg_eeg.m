% params = SO_params(1);
% cwd    = pwd;
% cd(params.SubjDir);
% 
% % combine_meg_eeg => detections
% % eeg events
% t       = te + (tg(1)-te(1));
% events  = [];
%   % The events structure can be used by visualizer to mark detections
% for k   = 1:length(detections)
%   events(k).label = eeg.sensor_info(k).label;
%   events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
%   events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
% end
% eeg_events = events; clear events
% % grad events
% load('SO_detections.mat','events');
% events(ng+1:ng+ne) = eeg_events;
% % save combined events structure
% save('SO_detections_grad_eeg.mat','events');

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


cd(cwd);