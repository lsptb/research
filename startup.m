% What about this? --GMF [4/22/14]
global BIOSIMROOT
BIOSIMROOT = pwd();
packages = {'/research/modeling', ...
'/research/analysis/mmil/matlab/timesurfer', ...
'/research/visualization/matlab/functions', ...
'/research/visualization/matlab/layouts', ...
'/research/analysis/mmil/matlab/master', ...
'/research/analysis/mmil/matlab/matext', ...
'/research/analysis/mmil/matlab/mmil_util', ...
'/research/analysis/mmil/matlab/csv', ...
'/research/analysis/mmil/matlab/fieldtrip-20080624_private', ...
'/research/analysis/mmil/matlab/eeglab_functions', ...
'/research/external/tallie', ...
'/research/external/tallie/funTA', ...
'/research/external/grant', ...
'/research/external/jsonlab', ...
'/research/external/circ_stats', ...
};

for p = 1:length(packages)
  addpath([BIOSIMROOT packages{p}]);
end

