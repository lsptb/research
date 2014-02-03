function [allspecs] = simstudy(spec,scope,variable,values,varargin)
% allspecs = simstudy(spec,scope,variable,values)
% allspecs = get_search_space(net,'E','multiplicity','[10:10:50]','sim_cluster','scc1.bu.edu');
% allspecs = get_search_space(net,{'{E,I}','E'},{'mechanisms','multiplicity'},{'{iNa,iK,ileak}','[10 20 30]'})

% Load spec
if ischar(spec) % this is a path to specification files
  if exist(spec,'dir')
    spec = loadspec(spec);
  else
    error('specification not find.');
  end
end

% parameters
spec.simulation = mmil_args2parms( varargin, ...
                   {  'logfid',1,[],...
                      'logfile',[],[],...
                      'sim_cluster_flag',1,[],...
                      'sim_cluster','scc2.bu.edu',[],...
                      'sim_qsubscript','qmatjobs_memlimit',[],...
                      'sim_driver','biosimdriver.m',[],...
                      'ProjName','model',[],...
                      'StudyName','study',[],...
                      'description',[],[],...
                      'dt',.01,[],...
                      'SOLVER','euler',[],...
                      'memlimit','8G',[],...
                      'batchdir',[],[],...
                   }, false);
                 
% get search space
spec.simulation.scope = scope;
spec.simulation.variable = variable;
spec.simulation.values = values;
allspecs = get_search_space(spec);

% log file
logfile=spec.simulation.logfile;
logfid=spec.simulation.logfid;
if ischar(logfile) && ~isempty(logfile), logfid = fopen(logfile,'w'); else logfid = 1; end

% system info
[o,host]=system('echo $HOSTNAME');
[o,home]=system('echo $HOME');
home=home(1:end-1); % remove new line character

if spec.simulation.cluster_flag % run on cluster
  % create batchdir
  if isempty(spec.simulation.batchdir)
    batchname = sprintf('%s_%s_%s',spec.simulation.ProjName,spec.simulation.StudyName,datestr(now,30)); 
    batchdir = sprintf('%s/batchdirs/%s',home,batchname);
    spec.simulation.batchdir=batchdir;
  else
    [fp,batchname]=fileparts(spec.simulation.batchdir);
    batchdir = spec.simulation.batchdir;
  end
  mkdir(batchdir);
  cwd = pwd; cd(batchdir);
  driverscript = which(spec.simulation.sim_driver);
  [fpath,scriptname,fext] = fileparts(driverscript);
  %allfiles = {driverscript spec.files{:}};
  %for i = 1:length(allfiles)
  %  unix(sprintf('cp %s .',allfiles{i}));
  %end
  % create jobs
  jobs={};
  for k=1:length(allspecs)
    modelspec=allspecs{k};
    specfile = sprintf('spec%g.mat',k);
    save(specfile,'modelspec');
    jobs{end+1} = sprintf('job%g.m',k);
    fileID = fopen(jobs{end},'wt');
    fprintf(fileID,'load(''%s'',''modelspec''); %s(modelspec);\n',specfile,scriptname);
    fprintf(fileID,'exit');
    fclose(fileID);
  end
  % create scriptlist.txt (list of jobs)
  fileID = fopen('scriptlist.txt', 'wt');
  for i = 1:length(jobs)
    [a,this] = fileparts(jobs{i});
    fprintf(fileID,'%s\n',this);
  end
  fclose(fileID);
  % submit the jobs
  cmd = sprintf('%s %s %s',spec.simulation.sim_qsubscript,batchname,spec.simulation.memlimit);
  fprintf(logfid,'executing: "%s" on cluster %s\n',cmd,spec.simulation.sim_cluster);
  if ~strmatch(host,spec.simulation.sim_cluster);
    % connect to cluster and submit jobs
    if 0
      [s,m] = unix(sprintf('ssh %s "%s"',spec.simulation.sim_cluster,cmd));
    end
  else
    % submit jobs on the current host
    [s,m] = unix(cmd);
  end
  % log errors
  if s, fprintf(logfid,'%s',m); end
  cd(cwd);
  fprintf(logfid,'Jobs submitted.\n');        
else
  % run on local machine
  for specnum = 1:length(allspecs) % loop over elements of search space
    modelspec = allspecs{specnum};
    fprintf(logfid,'processing simulation...');
    simdriver(modelspec);
    fprintf(logfid,'done (%g of %g)\n',specnum,length(allspecs));
  end
end

end % end main function
