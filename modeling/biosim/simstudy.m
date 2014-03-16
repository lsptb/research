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
                      'description',[],[],...
                      'dt',.01,[],...
                      'SOLVER','euler',[],...
                      'memlimit','8G',[],...
                      'batchdir',[],[],...
                      'rootdir',pwd,[],...
                      'override',[],[],...
                      'timelimits',[],[],...
                      'dsfact',[],[],...
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

% define output directory structure
scopes = cellfun(@(x)x.simulation.scope,allspecs,'uni',0);
vars = cellfun(@(x)x.simulation.variable,allspecs,'uni',0);
vals = cellfun(@(x)x.simulation.values,allspecs,'uni',0);
uniqscopes = unique(scopes);
outdirs={}; dirinds=zeros(size(allspecs));
for k=1:length(uniqscopes)
  scopeparts = regexp(uniqscopes{k},'[^\(\)]*','match');
  uniqscopeparts = unique(scopeparts);
  specind = find(strcmp(uniqscopes{k},scopes));
  varparts = regexp(vars{specind(1)},'[^\(\)]*','match');
  dirname = '';
  for j=1:length(uniqscopeparts)
    selind = find(strcmp(uniqscopeparts{j},scopeparts));
    dirname = [dirname '__' strrep(uniqscopeparts{j},',','+')];
    for i=1:length(selind)
      dirname = [dirname '-' strrep(varparts{selind(i)},'_','')];
    end
  end
  outdirs{end+1} = dirname(3:end);
  dirinds(specind)=k;
end
%   inds = strmatch(uniqscopes{k},scopes,'exact');
%   varparts = regexp(vars{inds(1)},'[^\(\)]*','match');
%   dirname = '';
%   for j=1:length(scopeparts)
%     dirname = [dirname '_' strrep(scopeparts{j},',','') '-' strrep(varparts{j},'_','')];
%   end
%   outdirs{end+1} = dirname(2:end);
%   dirinds(inds) = k;
% end
rootoutdir={}; prefix={};
timestamp = datestr(now,'yyyymmdd-HHMMSS');
for i=1:length(allspecs)
  rootoutdir{i} = fullfile(spec.simulation.rootdir,timestamp,outdirs{dirinds(i)});
  try
    scopeparts=regexp(allspecs{i}.simulation.scope,'[^\(\)]*','match');
    varparts = regexp(allspecs{i}.simulation.variable,'[^\(\)]*','match');
    valparts = regexp(allspecs{i}.simulation.values,'[^\(\)]*','match');
    uniqscopeparts = unique(scopeparts);
    pname='';
    for j=1:length(uniqscopeparts)
      selind = find(strcmp(uniqscopeparts{j},scopeparts));
      pname = [pname '__' strrep(uniqscopeparts{j},',','+')];
      for k=1:length(selind)
        pname = [pname '-' strrep(varparts{selind(k)},'_','') strrep(valparts{selind(k)},'.','pt')];
      end
    end
    prefix{i}=pname(3:end);
  catch
    tmp=regexp(allspecs{i}.simulation.description,'[^\d_].*','match');
    prefix{i}=strrep([tmp{:}],',','_');
  end
  fprintf(logfid,'%s: %s\n',rootoutdir{i},prefix{i});
end
% save allspecs(i) results in rootoutdir{i}

% system info
[o,host]=system('echo $HOSTNAME');
[o,home]=system('echo $HOME');
home=home(1:end-1); % remove new line character
cwd = pwd;

if spec.simulation.sim_cluster_flag % run on cluster
  % create batchdir
  if isempty(spec.simulation.batchdir)
    batchname = ['B' timestamp];
    batchdir = sprintf('%s/batchdirs/%s',home,batchname);
    spec.simulation.batchdir=batchdir;
  else
    [fp,batchname]=fileparts(spec.simulation.batchdir);
    batchdir = spec.simulation.batchdir;
  end
  mkdir(batchdir);
  cd(batchdir);
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
    fprintf(fileID,'load(''%s'',''modelspec''); %s(modelspec,''rootoutdir'',''%s'',''prefix'',''%s'',''cluster_flag'',1,''batchdir'',''%s'',''jobname'',''%s'');\n',specfile,scriptname,rootoutdir{k},prefix{k},batchdir,jobs{end});
    fprintf(fileID,'exit\n');
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
      [s,m] = system(sprintf('ssh %s "%s"',spec.simulation.sim_cluster,cmd));
    end
  else
    % submit jobs on the current host
    [s,m] = system(cmd);
  end
  % log errors
  if s, fprintf(logfid,'%s',m); end
  fprintf(logfid,'Jobs submitted.\n');        
else
  % run on local machine
  for specnum = 1:length(allspecs) % loop over elements of search space
    modelspec = allspecs{specnum};
    fprintf(logfid,'processing simulation...');
    biosimdriver(modelspec,'rootoutdir',rootoutdir{specnum},'prefix',prefix{specnum},'verbose',1,...
      'savefig_flag',0,'savedata_flag',0);
    fprintf(logfid,'done (%g of %g)\n',specnum,length(allspecs));
  end
end
cd(cwd);
end % end main function
