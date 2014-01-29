%  epoch_data = ntools_gen_cont_data(dat_file, experiment, ['channels', wh_channels_to_keep], ...
%       ['max_epochs', max_epochs], ['coords', coordinates])
%
%  Generate an epoch_data structure with 1 very long trial
%
%  Parameters:
%       dat_file: path to a .dat file
%       experiment: experiment structure
%
%  Optional keyword parameters:
%       channels: which channels to keep (default: all)
%       max_epochs: max num epochs (default: no limit)
%       coords: an nx3 array of electrode coordinates
%
%  Note:
%       Equivalent to calling ntools_gen_epoch_data(... 'is_continuous', 1);

function epoch_data = ntools_gen_cont_data(varargin)
    epoch_data = ntools_gen_epoch_data(varargin{:}, 'is_continuous', 1);
end