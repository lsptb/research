[o,r]=system('echo $HOME'); % get home directory
r=r(1:end-1); % remove new line
addpath(genpath(sprintf('%s/research/modeling',r)));
addpath(genpath(sprintf('%s/research/analysis/mmil/matlab/timesurfer',r)));
addpath(genpath(sprintf('%s/research/visualization/matlab/functions',r)));
addpath(genpath(sprintf('%s/research/visualization/matlab/layouts',r)));
addpath(sprintf('%s/research/analysis/mmil/matlab/master',r));
addpath(sprintf('%s/research/analysis/mmil/matlab/matext',r)); % eg. splitstr()
addpath(sprintf('%s/research/analysis/mmil/matlab/mmil_util',r));
addpath(sprintf('%s/research/analysis/mmil/matlab/csv',r));
addpath(sprintf('%s/research/analysis/mmil/matlab/fieldtrip-20080624_private',r));
addpath(sprintf('%s/research/analysis/mmil/matlab/eeglab_functions',r));
addpath(sprintf('%s/research/external/tallie',r));
addpath(sprintf('%s/research/external/tallie/funTA',r));
addpath(sprintf('%s/research/external/grant',r));
addpath(sprintf('%s/research/external/jsonlab',r));

%addpath /usr3/graduate/sherfey/research/analysis/mmil/matlab/mmil_util % eg. mmil_args2parms()
