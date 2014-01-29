function ts_write_session_log(infile)
%--------------------------------------------------------------------------
%  Purpose: To create a log of ts_session averaging outcomes
% 
%  Required input: proc_avg_data_1.mat file created by ts_session 
%                  **make sure to change directories to matfiles** 
%   
%  Created: 6/19/2009
%--------------------------------------------------------------------------  
 if(nargin<1)||(~ischar(infile));
    % either no input or input is not a string
    infile='session.mat'; %default filename
 elseif ~ischar(infile)
    % there is an input but it is not a string
    [filename pathname]=uigetfile({'*.mat;*.log;'},'Pick.');
    if isequal(filename,0) || isequal(pathname,0)
      return;
    end
    infile = fullfile(pathname,filename);     
end
if ~exist(infile,'file')
    fprintf('%s:File not found:%s\n',mfilename,infile);
    return;
end
% load sessionfile
load(infile); %session
%create file and write header
[fpath fname]=fileparts(infile);
%outfile = sprintf('%s/practice_session_log.txt',session(1).parms.rootoutdir);
% outfile = sprintf('%s /home/ayseo/Practice/%s_parms.log',fpath,fname);
outfile = fullfile(fpath,[fname '_parms.txt']);
fid = fopen(outfile,'at');
fprintf(fid,'Parameter log: %s\n',date);
fprintf(fid,'Session file: %s\n\n',infile);
fprintf(fid,'---------------------');
n = length(session);
for k=n;
    fun = session(n).name;
    id = session(n).function_id;
    fprintf(fid,'\n---------------------\nfunction:%s\n',fun);
    if ~isnan(id)
        fprintf(fid,'User assigned ID %g to this function call.\n',id);
    end
    fields = fieldnames(session(n).parms);
    for f=1:length(fields)
        thisparm = fields{f};
        thisval = session(n).parms.(thisparm);
        %str=''; %this will be for the string version of this val
        if ischar(thisval)
            str = thisval;
        elseif isnumeric(thisval)
            str = num2str(thisval);  
        elseif iscellstr(thisval)    %this may not be the correct function
            c  = cell2mat(thisval);
            str = mat2str(c);
        %elseif %continue checking for other datatypes
        end
        fprintf(fid, '%-30s: %s\n',thisparm,str);
    end %end loop over parameters
end %end loop over elements of the "session" structure array
%write the footer
fprintf('%s: parameter list compiled on %s\n',mfilename,date);
%end
fclose(fid);
fprintf('%s: saving file: %s\n',mfilename,outfile);
end

    