% % % [o,r]=system('echo $REPOS'); % get home directory
% % % r=r(1:end-1); % remove new line
% % % addpath(genpath(sprintf('%s/research/modeling',r)));
% % % addpath(genpath(sprintf('%s/research/analysis/mmil/matlab/timesurfer',r)));
% % % addpath(genpath(sprintf('%s/research/visualization/matlab/functions',r)));
% % % addpath(genpath(sprintf('%s/research/visualization/matlab/layouts',r)));
% % % addpath(sprintf('%s/research/analysis/mmil/matlab/master',r));
% % % addpath(sprintf('%s/research/analysis/mmil/matlab/matext',r)); % eg. splitstr()
% % % addpath(sprintf('%s/research/analysis/mmil/matlab/mmil_util',r));
% % % addpath(sprintf('%s/research/analysis/mmil/matlab/csv',r));
% % % addpath(sprintf('%s/research/analysis/mmil/matlab/fieldtrip-20080624_private',r));
% % % addpath(sprintf('%s/research/analysis/mmil/matlab/eeglab_functions',r));
% % % addpath(sprintf('%s/research/external/tallie',r));
% % % addpath(sprintf('%s/research/external/tallie/funTA',r));
% % % addpath(sprintf('%s/research/external/grant',r));
% % % addpath(sprintf('%s/research/external/jsonlab',r));
% % % addpath(sprintf('%s/research/external/circ_stats',r));
% % % 
% % % %addpath /usr3/graduate/sherfey/research/analysis/mmil/matlab/mmil_util % eg. mmil_args2parms()

% addpath([ pwd() '/modeling/')];
% addpath([ pwd() '/visualization/matlab/timesurfer']);
% addpath([ pwd() '/visualization/matlab/functions']);
% addpath([ pwd() '/visualization/matlab/layouts']);
% addpath([ pwd() '/analysis/mmil/matlab/master']);
% addpath([ pwd() '/analysis/mmil/matlab/matext']);
% addpath([ pwd() '/analysis/mmil/matlab/mmil_util']);
% addpath([ pwd() '/analysis/mmil/matlab/csv']);
% addpath([ pwd() '/analysis/mmil/matlab/fieldtrip-20080624_private']);
% addpath([ pwd() '/analysis/mmil/matlab/eeglab_functions']);
% addpath([ pwd() '/external/tallie']);
% addpath([ pwd() '/external/tallie/funTA']);
% addpath([ pwd() '/external/grant']);
% addpath([ pwd() '/external/jsonlab']);
% addpath([ pwd() '/external/circ_stats']);

% What about this? --GMF [4/22/14]
addpath(genpath(pwd));
