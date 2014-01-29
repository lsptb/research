function [data data_name] = ntools_elec_read_struct_data(s,name)
% read the data from 1*1 structure into a double array

k = 1;
for i = 1:length(name)
    for j = 1:length(s.(char(name{i})))
        data(k,:) = s.(char(name{i}))(j,:);
        data_name{k,:} = sprintf('%s%.2d',name{i},j);
        k = k+1;
    end
end