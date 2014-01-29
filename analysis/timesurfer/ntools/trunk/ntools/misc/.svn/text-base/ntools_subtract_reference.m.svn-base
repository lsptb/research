% epoch_data = ntools_subtract_reference(epoch_data, ref_channel_num)
%
% Subtracts out the data of a reference channel from all other channels,
% and returns the result as a new epoch_data structure.
%
% E.g.:
%   ref_channel_num = 28;
%   ref_epoch_data = ntools_subtract_reference(epoch_data, ref_channel_num);

function epoch_data = ntools_subtract_reference(epoch_data, ref_channel_num)
    num_conditions = length(epoch_data.epochs);
    
    for i=1:num_conditions
        epoch_data.epochs(i).data = bsxfun(@minus, epoch_data.epochs(i).data,  epoch_data.epochs(i).data(ref_channel_num,:,:));
    end
end