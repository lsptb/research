[o,r]=system('echo $HOME'); % get home directory
r=r(1:end-1); % remove new line
addpath(genpath(sprintf('%s/research/modeling',r)));
addpath(genpath(sprintf('%s/research/analysis/mmil/matlab/timesurfer',r)));
addpath(sprintf('%s/research/analysis/mmil/matlab/matext',r)); % eg. splitstr()
addpath(sprintf('%s/research/analysis/mmil/matlab/mmil_util',r));
addpath(sprintf('%s/research/analysis/mmil/matlab/fieldtrip-20080624_private',r));
addpath(sprintf('%s/research/analysis/mmil/matlab/csv',r));
%addpath /usr3/graduate/sherfey/research/analysis/mmil/matlab/mmil_util % eg. mmil_args2parms()
