function avg_data=btq_Ravgs(filenames,varargin)
%
% a function to read avg_data from filenames and return an avg data
% structure of the combined data.
% 
% optional argument is vector of events to pull from the files
%
tic
if ~iscellstr(filenames),
    filenames={filenames};
end

for f=1:length(filenames),
    dat=load(filenames{f});
    display('Removing .parms field from avg_data');
    if isfield(dat.avg_data,'parms')
        dat.avg_data=rmfield(dat.avg_data,'parms');
    end    
    if isempty(varargin),
        avg_data{f}=ts_data_selection(dat.avg_data);
    else
        avg_data{f}=ts_data_selection(dat.avg_data,'events',varargin{1});
    end
end

if f>1,
    avg_data=ts_combine_data(avg_data{:});
else
    avg_data=avg_data{1};
end


toc