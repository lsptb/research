function run_timesurfer(parameter_file,varargin)
% Required inputs:
%   parmfile or parms structure
%     parmfile - *.csv or *.xls file containing parameter specifications
%     parms - structure containing parameter specifications
%       (fields required for each function k): 
%         parms.function(k).name - name of the function or script
%         parms.function(k).input - name of the input timesurfer data structure (ex. avg_data)
%         parms.function(k).output - name of the output timesurfer data structure
%         parms.function(k).parms - function parameters
%       (optional fields):
%         parms.global - parameters to pass to all functions
%         parms.function(k).process_flags.cluster - [0|1] whether to run on the cluster
% Optional inputs:
%

% Created by Jason Sherfey on 21-Nov-2008
% Modified last by JSS on 08-Apr-2009
%   - replaced load_parameters_csvfile by load_parms
%   - added option to provide parm struct instead of parmameter_file
%   - changed data.opt to data.parms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STARTUP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
log_opt = [];

% load general user-specified parameters
if nargin==1 && ischar(parameter_file)
  % read csv or xls parmfile
  fprintf('%s: loading parameters\n',mfilename);
  if strfind(parameter_file,'.csv')
%     base = load_parameters_csvfile(parameter_file);
    base = load_parms(parameter_file);
  else
    base = load_parameters(parameter_file);
  end
elseif nargin==1 && isstruct(parameter_file)
  % use user-supplied parms structure
  base = parameter_file;
elseif nargin > 1 && ischar(parameter_file)
  % user-supplied function name and key/value pairs
  base.function.name  = parameter_file;
  base.function.parms = mmil_args2parms(varargin{:});
  [base.function.input base.function.output] = ts_function_info(parameter_file,'io');
end

% load subject parameters
if issubfield(base,'global.subjectfile') && ~isempty(base.global.subjectfile)
  subj = load_subject_parameters(base);
else
  subj = base;  % 1 subject
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1:length(subj)
  parms = subj(s);
  try 
    if parms.global.saveparms
      save(fullfile(parms.global.rootoutdir,'parms.mat'),'parms');
    end
  end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for f = 1:length(parms.function) % loop over functions
    %% FUNCTION INITIALIZATION
    fun            = parms.function(f).name;
    itype          = parms.function(f).itype;
    otype          = parms.function(f).otype;
    opt            = parms.function(f).parms;
    opt.function   = fun;
    opt.indata     = itype;
    opt.outdata    = otype;
    runfun_flag    = parms.function(f).process_flags.run_flag;
    script_flag    = parms.function(f).process_flags.script_flag;
    cluster_flag   = parms.function(f).process_flags.cluster_flag;
    load_flag      = parms.function(f).spec_flags.load_flag;
    save_flag      = parms.function(f).spec_flags.save_flag;    
%     fun           = parms.function(f).name;
%     itype         = parms.function(f).input;
%     otype         = parms.function(f).output;
%     opt           = parms.function(f).parms;
%     opt.function  = fun;
%     opt.indata    = itype;
%     opt.outdata   = otype;
%     runfun_flag   = parms.function(f).process_flags.run;
%     script_flag   = parms.function(f).process_flags.auto_script;
%     cluster_flag  = parms.function(f).process_flags.cluster;
    %hdrout_flag  = parms.function(f).process_flags.save_header;
    try loop = opt.loop_param; catch loop = {}; end
    [OPT,N] = loop_info(opt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i = 1:N % loop over multiple calls to the same function
      opt = OPT;
      try 
        opt.(loop) = opt.(loop){i};
        try opt.statistics = opt.statistics{i}; end    
      end
      if isfield(opt,loop) && length(opt.(loop))==2
        fileid = sprintf('_%g_%g',opt.(loop)(1),opt.(loop)(2));
      elseif isfield(opt,loop) && length(opt.(loop))==1
        fileid = sprintf('_%g',opt.(loop)(1));
      else
        fileid = '';
      end    
      load_flag = 1;
      save_flag = 1;
      if isfield(opt,'hdr') && ~isempty(opt.hdr)
        load_flag = 0;
        itype = 'hdr';
        if ~isfield(opt,'datafile') || isempty(opt.datafile)
          opt.datafile = opt.hdr;
        end
      end
%% i=0 o=0 (wrapper)          
      if ~any(itype)
        load_flag = 0;        
        if ~any(otype)      
          save_flag = 0;
          fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s\n',...
            mfilename,s,length(subj),f,length(parms.function),i,N,fun);  
          [opt,instr] = wrapper_check(opt,fun,load_flag,save_flag);
          args = mmil_parms2args(opt);
          cmdstr = sprintf('%s(%s);',fun,instr);
          if cluster_flag
            if i==1, funcall = {}; batchdir = sprintf('%s_%s',fun,datestr(now,30)); end
            [funcall{i},outpath] = write_script(1,opt,fun,cmdstr,fileid,batchdir); 
          end    
          if script_flag,  write_script(2,opt,fun,cmdstr,fileid);   end
          if runfun_flag && ~cluster_flag
              if any(findstr(cmdstr,'parms')),cmdstr=strrep(cmdstr,'parms','opt'); end
              eval(cmdstr); 
          end
%% i=0 o=1 (fun loads data)          
        else
          fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s = %s(args{:})\n',...
            mfilename,s,length(subj),f,length(parms.function),i,N,otype,fun);              
          cmdstr = sprintf('%s = %s(args{:});',otype,fun);
        end
      else
%% i=1 o=0 (fun saves results)        
        if ~any(otype)
          save_flag = 0;
          cmdstr = sprintf('%s(%s,args{:});',fun,itype);
%% i=1 o=0 (i=raw)                    
          if strcmp(itype,'raw')  
            % for future use
          end
        else    
%% i=1 o=1 (i=raw)                
          if strcmp(itype,'raw')
            % for future use
%% i=1 o=1 (i=datafile) (fun loads data)            
          elseif strcmp(itype,'datafile')
            load_flag = 0;
            opt.indata = 0;
            fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s = %s(args{:})\n',...
              mfilename,s,length(subj),f,length(parms.function),i,N,otype,fun);              
            cmdstr = sprintf('%s = %s(args{:});',otype,fun);
%% i=1 o=1 (i=hdr) (fun loads data) not figure      
          elseif strcmp(itype,'hdr') && ~strcmp(otype,'figure')
            fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s = %s(hdr,args{:})\n',...
              mfilename,s,length(subj),f,length(parms.function),i,N,otype,fun);              
            cmdstr = sprintf('%s = %s(getfield(load(opt.hdr),''hdr''),args{:});',otype,fun);            
%% i=1 o=1 (i=hdr) (fun loads data) and figure      
          elseif strcmp(itype,'hdr') && strcmp(otype,'figure')
            fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s = %s(hdr,args{:})\n',...
              mfilename,s,length(subj),f,length(parms.function),i,N,otype,fun);              
            cmdstr = sprintf('%s(getfield(load(opt.hdr),''hdr''),args{:});',fun);                        
%% i=1 o=1 (o=figure) (plotting function)             
          elseif strcmp(otype,'figure')         
            save_flag = 2; cluster_flag = 0;
            fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s(%s)\n',...
              mfilename,s,length(subj),f,length(parms.function),i,N,fun,itype);              
            cmdstr = sprintf('%s(%s,args{:});',fun,itype);
%% i=1 o=1 (timesurfer function)                 
          else       
            fprintf('%s: subject %g of %g, function %g of %g, iteration %g of %g: %s = %s(%s)\n',...
              mfilename,s,length(subj),f,length(parms.function),i,N,otype,fun,itype);              
            cmdstr = sprintf('%s = %s(%s,args{:});',otype,fun,itype);  
          end
        end
      end
      if load_flag || save_flag
        % make sure datafiles and output path exist
        [opt,err] = check_files_and_paths(opt,load_flag,save_flag);
        if err==2, continue; end     
        % write scripts and run function
        if cluster_flag
          if ~exist('batchdir','var'), funcall = {}; batchdir = sprintf('%s_%s',fun,datestr(now,30)); end
          [funcall{i},outpath] = write_script(3,opt,fun,cmdstr,fileid,batchdir); 
        end
        if script_flag,  write_script(4,opt,fun,cmdstr,fileid);   end
        if runfun_flag && ~cluster_flag
          if load_flag
            % load data
            datafile = opt.datafile;
            comopt 	 = []; 
            for i = 1:length(datafile)
              fprintf('%s: loading data file %g of %g: %s\n',mfilename,i,length(datafile),datafile{i});
              tmpdat = getfield(load(datafile{i},itype),itype); 
              % prevent dissimilar structures error from inconsistent opts
              if isfield(tmpdat,'opt'), tmpdat.parms = tmpdat.opt; tmpdat = rmfield(tmpdat,'opt'); end
              if ~isfield(tmpdat,'parms'), tmpdat.parms = []; end
              if isempty(comopt) && ~isempty(tmpdat.parms)
                comopt = tmpdat.parms;
              end
              if i > 1, eval(sprintf('tmpdat = orderfields(tmpdat,%s(1));',itype)); end
              eval(sprintf('%s(i) = tmpdat;',itype));
              clear tmpdat;
              eval(sprintf('try hdr(i) = ts_get_header(%s(i)); end',itype));
            end
            eval(sprintf('[%s.parms] = deal(comopt);',itype));
            eval(sprintf('if length(%s)>1, %s=ts_combine_data(%s); end',itype,itype,itype));
            try 
              if parms.global.saveheader && exist('hdr','var')
    %                 save(fullfile(parms.global.rootoutdir,sprintf('%s_inhdr.mat',fun)),'hdr');
                save(fullfile(opt.rootoutdir,sprintf('%s_inhdr.mat',fun)),'hdr');
              end
            end           
            tmpargs = mmil_parms2args(opt);
            eval(sprintf('[%s opt] = ts_data_selection(%s,''opt'',opt,tmpargs{:});',itype,itype));
          end
          % call the timesurfer function
          fprintf('%s: calling the timesurfer function\n',mfilename);
          args = mmil_parms2args(opt);
          eval(cmdstr);
          % save the results
          if save_flag == 1
            eval(sprintf('try opt.previous = %s.parms; end',itype));
            eval(sprintf('%s.parms = opt;',otype));
%             cmdstr = sprintf('%s = outdata;', otype);
%             eval(cmdstr);
            try eval(sprintf('save(%s.parms.filename{1},''%s'');',otype,otype));
            catch
              keyboard
            end
            % save hdr ==> get hdr from otype; and save as ~ filename{1}.hdr.mat
          elseif save_flag ==2
              save_figure(args{:});
          end
        end
%         if isfield(log_opt,fun)
%           log_opt(s).(fun)(end+1) = opt;
%         else
%           log_opt(s).(fun) = opt;
%         end
        % free up some memory
        if ischar(itype), eval(sprintf('clear %s',itype)); end
        if ischar(otype), eval(sprintf('clear %s',otype)); end 
      end
    end % loop over multiple calls to the same function
    if cluster_flag
      write_scriptlist(outpath,funcall);
%       fprintf('---------------------------------------------------------------------------------------\n');
%       fprintf('type at prompt: qmatjobs %s\nthen type exit to resume processing\n',batchdir)
%       fprintf('---------------------------------------------------------------------------------------\n');
      if isfield(opt,'cluster') && ischar(opt.cluster)
        mmilcluster = opt.cluster;
      else
        switch mod(f,3)
          case 0
            mmilcluster = 'mmilcluster';
          case 1
            mmilcluster = 'mmilcluster2';
          case 2
            mmilcluster = 'mmilcluster3';
        end
      end
%       eval(sprintf('!ssh -XY %s',mmilcluster));      
      fprintf('submitting job to cluster ''%s'':\n%s\n',mmilcluster,batchdir);
      jobstring = sprintf('qmatjobs %s',batchdir);
      eval(sprintf('!ssh %s %s',mmilcluster,jobstring));
    end % end-if cluster
    toc
  end   % end-loop over different functions
  toc
end     % end-loop over subjects

% batch processing of multiple subjects if auto_script and no individual
% cluster nor local processing (make subject_batches_flag at beginning)

% for each subjectdir
%   write_scriptlist(subjectdir/scripts,funcall)
%   copy scriptlist and all funcall from subjectdir/scripts to
%   user/batchdir/subjectid_uid
%   ssh to 1 of 3 clusters and prompt user
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEANUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save reports / logs
    
    if ~isempty(log_opt)
      % print text file with info in log_opt: session summary
      writesummary(log_opt);
      % save log_opt
      save('session_parameters.mat','log_opt');
    end

% remove any temporary files

% memory dump
clear all
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [opt,err] = check_files_and_paths(opt,load_flag,save_flag)
% err: 0-ok, 1-error, 2-skip function

[opt,err,str] = check_files(opt,load_flag,save_flag);
if ~err && save_flag, [err,str] = make_paths(opt); end

if err==1, error('%s: %s',mfilename,str);   end
if err==2, fprintf('%s: %s',mfilename,str); end

function [opt,err,str] = check_files(opt,load_flag,save_flag)
% check input and output files
err = 0; str = ''; 
if ~isfield(opt,'overwrite'), opt.overwrite = 0; end
if ~isfield(opt,'prompt'), opt.prompt = 0; end
% if load_flag
  % timesurfer data input
  if ~isfield(opt,'datafile') || isempty(opt.datafile)
    opt = find_datafiles(opt);
%     opt = find_timesurfer_files(opt);
    if isempty(opt.datafile)
      str = sprintf('could not find data files\nSkipping %s\n',opt.function);
      err = 2; % skip function
      return;
    end  
  else
    % datapath
    if ~isfield(opt,'datapath')
      try   opt.datapath = opt.rootindir;
%       try 
%         [opt.datapath jnk] = fileparts(opt.datafile);
      catch
        opt.datapath = pwd;
      end
    end  
  end
  if ~iscell(opt.datafile), opt.datafile = {opt.datafile}; end

  if isempty(opt.datafile)
    str = sprintf('could not find %s for these events\nskipping %s\n',opt.indata,opt.function);
    err = 2;
    return;
  end

  % verify that datafiles exist
  for i = 1:length(opt.datafile)
    file = opt.datafile{i};
    tmp = file;
    if ~exist(file,'file')
      try tmp = fullfile(opt.datapath,file); end
    end
    if ~exist(tmp,'file')
      try tmp = fullfile(opt.datapath,opt.inpath,file); end
    end
    if ~exist(tmp,'file')
      try tmp = fullfile(opt.rootindir,file); end
    end
    if ~exist(tmp,'file')
      try tmp = fullfile(opt.rootindir,opt.inpath,file); end
    end
    opt.datafile{i} = tmp;
    clear tmp;
  end

  % make list of other input files & verify that they exist
  fprintf('%s: verifying files exist...',mfilename);
  filetype = {'badchanfile','layoutfile','rejectfile','baselinefile','statfile',...
              'badchan_file','layout_file','reject_file','baseline_file','stat_file'};
  f      = fieldnames(opt);
  ind    = find(isfield(opt,filetype));
  tmp_name = {};
  for j = 1:length(ind)
    try tmp_name_ = opt.(filetype{ind(j)}); catch tmp_name_ = []; end
    if isempty(tmp_name_), continue; end
    charflag = 0;
    if ~iscell(tmp_name_), tmp_name_ = {tmp_name_}; charflag = 1; end
    for k = 1:length(tmp_name_)
      tmp = tmp_name_{k};
      if ~exist(tmp,'file')
        try tmp = fullfile(opt.datapath,tmp_name_{k}); end
      end
      if ~exist(tmp,'file')
        try tmp = fullfile(opt.datapath,opt.inpath,tmp_name_{k}); end
      end		
      if ~exist(tmp,'file')
        try tmp = fullfile(opt.rootindir,tmp_name_{k}); end;
      end     
      if ~exist(tmp,'file')
        try tmp = fullfile(opt.rootindir,opt.inpath,tmp_name_{k}); end;
      end
      if ~exist(tmp,'file') && opt.prompt
        fprintf('%s does not exist: %s\n',filetype{ind(j)},tmp);
        title = sprintf('%s not found for %s.  select a file.\n',filetype{ind(j)},opt.function);
        [opt.datafile,opt.datapath] = uigetfile('*.mat',title);    
        tmp = fullfile(pathname,filename);    
      end
      if ~exist(tmp,'file')
        err = 2; % skip function
        str = sprintf('\nfailed to find input file: %s\nskipping %s\n',tmp,opt.function);
        return;
      end    
      tmp_name_{k} = tmp;
      if charflag
        opt.(filetype{ind(j)}) = tmp;
      else
        opt.(filetype{ind(j)}){k} = tmp{:};
      end
    end
    sz = length(tmp_name);  
    tmp_name(sz+1:sz+length(tmp_name_)) = tmp_name_;
  end
  opt.files = tmp_name;
  if isempty(opt.files), opt.files = []; end
  fprintf('done\n');
% end

if save_flag
  % get output filenames
  if ~issubfield('opt','previous.filename') && ~isempty(opt.datafile)
    opt.previous.filename = opt.datafile{1}; 
  else
    opt = find_datafiles(opt);
    opt.previous.filename = opt.datafile{1};
  end
  opt = ts_make_output_filename(opt);
  if ~iscell(opt.filename), opt.filename = {opt.filename}; end

  % do output files already exist.  if so, check opt.overwrite
  for i = 1:length(opt.filename)
      if exist(opt.filename{i},'file')
          if opt.overwrite, 
              warning('%s: file will be overwritten: %s\n',mfilename,opt.filename{i}); 
          else
              err = 2;
              str = sprintf('output file already exists: %s\nSkipping %s\n',opt.filename{i},opt.function);
          end
      end
  end
end

function [err,str] = make_paths(opt)
% makes output paths if they don't already exist
err = 0; str = '';
if ~iscell(opt.filename), opt.filename = {opt.filename}; end
[pathstr,name,ext,versn] = fileparts(opt.filename{1});
if ~exist(pathstr,'dir'),
    fprintf('making output directory: %s\n',pathstr);
    unix(['mkdir -p ' pathstr]);
end
opt.savepath = pathstr;

function [opt,N] = loop_info(opt)
if ~isfield(opt,'loop_param') || isempty(opt.loop_param), opt.loop_param = ''; N=1; return; end
p = opt.(opt.loop_param);
if ~iscell(p), p = {p}; end
N = length(p);
  
function write_scriptlist(outpath,funcall)
outfile = fullfile(outpath,'scriptlist.txt');
fprintf('%s: writing script list for cluster computing: %s\n',mfilename,outfile);
fid = fopen(outfile,'wt');
for i = 1:length(funcall), fprintf(fid,'%s\n',funcall{i}); end
fclose(fid);

function varargout = write_script(type,opt,fun,cmdstr,fileid, varargin)
if ~isfield(opt,'filename')
	opt.filename{1} = sprintf('run_%s.mat',fun);
end
switch type
  case 1      % cluster computing script: timesurfer wrapper
    batchdir = varargin{1};
    [e u] = unix('whoami');
    outpath = sprintf('/home/%s/batchdirs/%s',u(1:end-1),batchdir);
    funcall = sprintf('run_%s%s',fun,fileid);
    outfile = sprintf('%s/%s.m',outpath,funcall); 
    auto_script('options',opt,'function',fun,'cmdstr',cmdstr,'outfile',outfile,'fileid',fileid);   
    varargout{1} = funcall;
    varargout{2} = outpath;
  case 2      % not for cluster computing: timesurfer wrapper
    auto_script('options',opt,'function',fun,'cmdstr',cmdstr,'fileid',fileid);
  case 3      % cluster computing script: timesurfer function    
    batchdir = varargin{1};
    [e u] = unix('whoami');
    outpath = sprintf('/home/%s/batchdirs/%s',u(1:end-1),batchdir);
    funcall = sprintf('run_%s%s',fun,fileid);
    outfile = sprintf('%s/%s.m',outpath,funcall); 
    auto_script('options',opt,'function',fun,'cmdstr',cmdstr,'outfile',outfile,...
      'fileid',fileid,'filename',opt.filename{1});   
    varargout{1} = funcall;
    varargout{2} = outpath;
  case 4      % not for cluster computing: timesurfer function
    auto_script('options',opt,'function',fun,'cmdstr',cmdstr,'fileid',fileid,'filename',opt.filename{1});
end
function save_figure(varargin)
parms = mmil_args2parms(varargin,...
						{'format','eps',[],...
						 'rootoutdir',pwd,[],...
						 'savepath',[],[],...
                         'filename',[],[],...
						 'save',1,{1,0},...
						 'close',0,{1,0},...
                         'prefix','proc',[],...
                         'overwrite',0,{1,0}...
						},false);

if isempty(parms.filename)
    parms.filename{1} = sprintf('%s/images/%s_%s',parms.rootoutdir,parms.prefix,datestr(now,30));
else
    [pathstr,fname] = fileparts(parms.filename{1});
    parms.filename{1} = sprintf('%s/images/%s',parms.rootoutdir,fname);
end
[pathstr,fname,ext,versn] = fileparts(parms.filename{1});
if ~exist(pathstr,'dir'),
    fprintf('making output directory: %s\n',pathstr);
    unix(['mkdir -p ' pathstr]);
end
if strcmp(parms.format,'jpg'), parms.printcmd = {'-djpeg'}; end
if strcmp(parms.format,'eps'), parms.printcmd = {'-depsc','-tiff','-r150'}; end
if strcmp(parms.format,'tif'), parms.printcmd = {'-dtiff'}; end
fprintf('%s: saving figure: %s\n',mfilename,parms.filename{1});
if parms.save, print(gcf,parms.printcmd{:},parms.filename{1}); end
if parms.close, close; end

function writesummary(log_opt)
fprintf('%s: writing session summary and saving session parameters.\n',mfilename);
fid = fopen('session_parameters.log','wt');
fprintf(fid,'parameter record for successful functions (cells and structs are not shown)\n');
for i = 1:length(log_opt)         % subjects
  subj = log_opt(i);
  funs = fieldnames(subj);
  for j = 1:length(funs)          % functions
    fun = funs{i};
    funparm = subj.(fun);
    for k = 1:length(funparm)     % function iterations
      fprintf(fid,'\n---------------------------------------------------------------------------\n');
      fprintf(fid,'subject %d, iteration %d of function ''%s''\n',i,k,fun);
      fprintf(fid,'---------------------------------------------------------------------------\n\n');
      opt   = funparm(k);
      parms = fieldnames(opt);
      for p = 1:length(parms)
        parm = parms{p};
				pclass = class(opt.(parm));
        switch pclass
          case 'logical'
            value = opt.(parm);
            cchar = '%d';
          case 'char'
            value = opt.(parm);
            cchar = '%s';
          case 'struct'
            value = pclass;
            cchar = '%s';
					case 'cell'
						if length(opt.(parm))==1 && ischar(opt.(parm){1})
							value = opt.(parm){1};
						else
							value = pclass;
						end
						cchar = '%s';
          otherwise
            if isnumeric(opt.(parm))
              value = num2str(opt.(parm));
              cchar = '[%s]';
            else
              continue;
            end
        end
        cmdstring = sprintf('%s = %s;\\n',parm,cchar);
        eval(sprintf('fprintf(fid,''%s'',value);',cmdstring));
      end
    end
  end
end    
fclose(fid);

function [opt,instr] = wrapper_check(opt,fun,load_flag,save_flag)
if strcmp(fun,'ts_process_fif_data')
  [opt,err] = check_files_and_paths(opt,load_flag,save_flag);
  instr = ['parms.datafile,args{:}'];
elseif strcmp(fun,'ts_iEEG_ProcessEEG')
  [opt,err] = check_files_and_paths(opt,load_flag,save_flag);
%   args = mmil_parms2args(opt);
  instr = 'args{:}';
else
  instr = 'args{:}';
end
