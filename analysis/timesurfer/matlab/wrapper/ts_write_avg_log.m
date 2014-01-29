function ts_write_average_log(infile)
%--------------------------------------------------------------------------
%  Purpose: To create a log of ts_session averaging outcomes
% 
%  Required input: proc_avg_data_1.mat file created by ts_session 
%                  **make sure to change directories to matfiles** 
%   
%  Created: 6/19/2009
%--------------------------------------------------------------------------  
if(nargin<1)||(~ischar(infile));
  [filename pathname]=uigetfile({'*.mat;*.log;'},'Pick.');
  if isequal(filename,0) || isequal(pathname,0)
    return;
  end   
  infile = fullfile(pathname,filename);  
end

if ~exist(infile,'file')
    fprintf('%s: File not found:%s\n',mfilename,infile);
    return;
end
% load sessionfile
load(infile); %avg_data
%create file and write header
[fpath fname]=fileparts(infile);
%outfile = sprintf('%s/practice_session_log.txt',session(1).parms.rootoutdir);
% outfile = sprintf('%s /home/ayseo/Practice/%s_parms.log',fpath,fname);
outfile = fullfile(fpath,[fname '_log.txt']);
fid = fopen(outfile,'wt');
fprintf(fid,'Averaging parameters log: %s\n',date);
fprintf(fid,'Average file: %s\n\n',infile);
fprintf(fid,'-----------------------\n');

fprintf(fid,'number of sensors:  %g\n',avg_data.num_sensors);
fprintf(fid,'sampling frequency: %g\n',avg_data.sfreq);

% number of events used:
numevnts = length([avg_data.averages.event_code]);
      
for k=1:numevnts
    num_trials = avg_data.averages(k).num_trials;
     
    if ischar(num_trials);
        trials = num_trials;
    elseif isnumeric(num_trials);
        trials = num2str(num_trials);
    elseif iscellstr(num_trials);
        mat    = cell2mat(num_trials);
        trials = mat2str(c);
    elseif isstruct(thisval);
        cell   = struct2cell(num_trials);
        mat    = cell2mat(cell);
        trials = mat2str(mat);
    end
    event_codes = avg_data.averages(k).event_code;
    events      = num2str(event_codes); 
    fprintf(fid,'\n-----------------------\n');
    fprintf(fid,'\nevent number: %-20s\n\n',events);
    fprintf(fid, 'number of trials: %-20s\n\n', trials);

    num_rejects = avg_data.averages(k).num_rejects;
    rej_info    = fieldnames(num_rejects);
    fprintf(fid,'number of rejects:\n sensor type:         number rejected:\n');
    for r =1:length(rej_info);
      rej_fld = (rej_info{r});
      rej_val = num_rejects.(rej_fld);
      if ischar(rej_val);
        rejects = rej_val;
        elseif isnumeric(rej_val);
            rejects = num2str(rej_val);
        elseif iscellstr(rej_val);
            mat = cell2mat(rej_val);
            rejects = mat2str(c);
        elseif isstruct(thisval);
        cell = struct2cell(rej_val);
        mat  = cell2mat(cell);
        rejects = mat2str(mat);
      end
      fprintf(fid,' %-20s %-20s\n',rej_fld,rejects);
    end    
end

fclose(fid);
fprintf('%s: saving file: %s\n',mfilename,outfile);
end
