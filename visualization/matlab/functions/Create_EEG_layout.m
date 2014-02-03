% Create EEG layout file

xx = [.407 .5 .593 .325 .405 .5 .595 .675 .25 .31 .37 .43 .5 .57 .63 .69 .75 .19 ...
      .26 .325 .38 .44 .5 .56 .62 .675 .74 .81 .18 .245 .31 .38 .44 .5 .56 ...
      .62 .69 .755 .82 .185 .25 .315 .38 .44 .5 .56 .62 .685 .75 .815 .195 ...
      .26 .325 .38 .44 .5 .56 .62 .675 .74 .805 .31 .41 .5 .59 .69 .41 .5 .59 .5];
 
yy = [.89 .89 .89 .85 .81 .78 .81 .85 .78 .76 .74 .72 .71 .72 .74 .76 .78 .67 ...
      .66 .65 .64 .63 .64 .63 .64 .65 .66 .67 .55 .55 .55 .55 .55 .55 .55 ...
      .55 .55 .55 .55 .43 .44 .45 .46 .47 .46 .47 .46 .45 .44 .43 .3 .32 ...
      .35 .37 .39 .39 .39 .37 .35 .32 .3 .25 .29 .31 .29 .25 .21 .23 .21 .15];  
 
tmp       = ts_data_selection(data,'chantype','eeg','badlabels',{'EOG 061','EOG 062','ECG 063','EMG 064'});

eeglabels = {tmp.sensor_info.label};
% eeglabels = ...
%   {'FP1' 'FPz' 'FP2' 'AF7' 'AF3' 'AFz' 'AF4' 'AF8' 'F7' 'F5' 'F3' 'F1' 'Fz' ...
%    'F2' 'F4' 'F6' 'F8' 'FT9' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'FC6' ...
%    'FT8' 'FT10' 'T9' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'T10' 'TP9' ...
%    'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P9' 'P7' 'P5' ...
%    'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'EOG061' 'EOG062' 'ECG063' 'EMG064' ...
%    'P10' 'PO7' 'PO3' 'POz' 'PO4' 'PO8' 'O1' 'Oz' 'O2' 'Iz'};
% ii = [1:60 65:74]; eeglabels = eeglabels(ii);

lognum    = [tmp.sensor_info.lognum];

for k = 1:length(xx)
  fprintf('%g\t%g\t%g\t%g\t%g\t%s\n',lognum(k),xx(k),yy(k),.04,.04,eeglabels{k});
end

