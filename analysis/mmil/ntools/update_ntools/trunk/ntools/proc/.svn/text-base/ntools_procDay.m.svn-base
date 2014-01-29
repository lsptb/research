function ntools_procDay(day,startOn)
%
%  ntools_procDay(day, startOn)
%
%   Inputs:  DAY = String.  Recording day.  For example '080926';
%
% StartOn specifies which step of procDay to begin with (can be used to
% complete a procDay that was interrupted by an error.)  See
% procErrorLog.txt in the day folder for error description.  Steps are :
% 
% 1.procDat             
% 2.procSp              
% 3.procMlfp           
% 4.procDisplay        
% 5.procDayEvents       
% 6.procState           
% 7.saveTrials          
% 8.procLowPassEye	  
% 9.procCalibrateHand   
% 10.procCalibrateJoy   
% 11.saveTrials         
% 12.saveStates         
% 13.procChase          
% 14.saveChase           
global EXPERIMENTDIR MONKEYRECDIR
if(~isdir([EXPERIMENTDIR '/' day]))
    try
        disp('Copying day folder to raid drive.  This might take a while.');
        unix(['cp -rv ' EXPERIMENTRECDIR '/' day ' ' EXPERIMENTDIR '/'],'-echo');
        unix(['chown bijanadmin:bijanadmin ' EXPERIMENTDIR '/' day]);
        disp([EXPERIMENTDIR '/' day ' ownership set'])
    catch
        disp('failed to copy recording folder:')
        disp(lasterror.message);
        return;
    end
else
    disp('Directory already copied to raid drive')
end
olddir = pwd;
recs = dayrecs(day);
% nRecs = length(recs);

% if nargin < 2
%     num = [1,nRecs];
% elseif ischar(rec) 
%     num = [find(strcmp(recs,rec)),find(strcmp(recs,rec))];
% elseif length(rec)==1 
%     num = [rec,rec];
% elseif length(rec)==2
%     num = rec;
% end

if nargin < 2 || startOn > 14 || startOn < 1
    startOn = 1;
end

% load([EXPERIMENTDIR '/' day '/' recs{1} '/rec' recs{1} '.Rec.mat']);
% if nargin < 2; ch{1} = 1:Rec.Ch(1); ch{2} = 1:Rec.Ch(2); end
stuffToDo = cell(0);

stuffToDo{end+1} = 'procDat'; %1
stuffToDo{end+1} = 'procSp';  %2
stuffToDo{end+1} = 'procMlfp';%3
stuffToDo{end+1} = 'procDisplay'; %4
stuffToDo{end+1} = 'procDayEvents'; %5
stuffToDo{end+1} = 'procState';%6
stuffToDo{end+1} = 'saveTrials'; %7
stuffToDo{end+1} = 'procLowPassEye'; %8
stuffToDo{end+1} = 'procCalibrateHand'; %9
stuffToDo{end+1} = 'procCalibrateJoy'; %10
stuffToDo{end+1} = 'saveTrials'; %11
stuffToDo{end+1} = 'saveStates'; %12
stuffToDo{end+1} = 'procChase'; %13
stuffToDo{end+1} = 'saveChase';%14

for i=startOn:length(stuffToDo)
    try
        feval(stuffToDo{i},day);
    catch
        log = 'ERROR!\nSucceeded on the following steps:';
        if(i>1)
            for j=1:i-1;
                log = [log  '\n' stuffToDo{j}];
            end
        end
        err = lasterror;
        log = [log '\nFailed on step ' stuffToDo{i} ' for reason:\n' err.message];
    
        s = err.stack;
        for j=1:length(s);
            log = [log '\n ...in ' s(j).name];
            log = [log '\n' s(j).file];
            log = [log '\n...at line ' num2str(s(j).line)];
        end 
     
        
        fid = fopen([EXPERIMENTDIR '/' day '/procErrorLog.txt'],'w');
        msg = sprintf(log);
        fwrite(fid,msg);
        fclose(fid);
        disp(msg);
        disp('Log saved.  Aborting.');
        cd(olddir);
        return;
    end
    cd(olddir);
end



% delFiles(day);  In the past, when I automatically deleted files I deleted entire recording days
%                   by accident.  I think delFiles should be done manually unless there is
%                   a procedure for automatically determining that the preprocessing was successful.
