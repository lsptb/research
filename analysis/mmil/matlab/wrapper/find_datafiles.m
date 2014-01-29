function opt = find_timesurfer_files(opt)
% datapath
if ~isfield(opt,'prompt'), prompt=0; else prompt=opt.prompt; end
if ~isfield(opt,'rootoutdir') || isempty(opt.rootoutdir), opt.rootoutdir = pwd; end
if ~isfield(opt,'rootindir') || isempty(opt.rootindir), opt.rootindir = opt.rootoutdir; end
if ~isfield(opt,'datapath') || isempty(opt.datapath)
  opt.datapath = opt.rootindir;
  if isfield(opt,'inpath') && ~isempty(opt.inpath)
    opt.datapath = fullfile(opt.datapath,opt.inpath);
  end
end
if (~isfield(opt,'indata')  || isempty(opt.indata))  && (isfield(opt,'itype') && ~isempty(opt.itype))
  opt.indata = opt.itype;
end

% backward compatibility
if isfield(opt,'datasubstring'), opt.datastring = opt.datasubstring; end
if isfield(opt,'datatype'),      opt.datastruct = opt.datatype;      end
if isfield(opt,'dataexclude'),   opt.dataskip = opt.dataexclude;     end
if isfield(opt,'dataignore'),    opt.dataskip = opt.dataignore;      end

% given: rootoutdir, rootindir, datapath
datafile = {};

fprintf('%s: searching for datafiles in %s\n',mfilename,opt.datapath);  
if any(isfield(opt,{'dataprefix','datasuffix','datastring','dataskip','datastruct'}))
  files = dir(opt.datapath);
  files = {files.name};
  k = 1;
  for f = 1:length(files)
    [pathstr,name,ext,versn] = fileparts(files{f}); 

    if isfield(opt,'ext') && ~isempty(opt.ext) 
      if strcmp(ext,['.' opt.ext])
        datafile{k} = files{f}; 
      else
        if ~isempty(datafile) && ischar(datafile{end}), datafile{k} = {}; end
        continue;
      end
    end
    if isfield(opt,'dataprefix') && ~isempty(opt.dataprefix)
      if strmatch(opt.dataprefix,name)
        datafile{k} = files{f};
      else
        if ~isempty(datafile) && ischar(datafile{end}), datafile{k} = {}; end
        continue;
      end
    end
    if isfield(opt,'datasuffix') && ~isempty(opt.datasuffix)
      if strmatch(fliplr(opt.datasuffix),fliplr(name))
        datafile{k} = files{f};
      else
        if ~isempty(datafile) && ischar(datafile{end}), datafile{k} = {}; end
        continue;
      end
    end
    if isfield(opt,'datastring') && ~isempty(opt.datastring)
      if any(strfind(name,opt.datastring))
        datafile{k} = files{f};
      else
        if ~isempty(datafile) && ischar(datafile{end}), datafile{k} = {}; end
        continue;
      end
    end
    if isfield(opt,'dataskip') && ~isempty(opt.dataskip)
      if any(strfind(name,opt.dataskip))
        if ~isempty(datafile) && ischar(datafile{end}), datafile{k} = {}; end
        continue;
      else
        datafile{k} = files{f};
      end
    end
    if isfield(opt,'datastruct') && any(opt.datastruct)
      if strcmp(ext,'.mat') && any(strcmp(who('-file',fullfile(opt.datapath,files{f})),opt.datastruct))
        datafile{k} = files{f};
      else
        if ~isempty(datafile) && ischar(datafile{end}), datafile{k} = {}; end
        continue;
      end
    end
    k = k + 1;
  end
  if ~isempty(datafile) && isempty(datafile{end}) && length(datafile) > 1
    datafile = datafile(1:end-1);
  end
else
  datafile = {};
  if isfield(opt,'ext') && ~isempty(opt.ext)
    fprintf('%s: searching for %s datafile in %s\n',mfilename,opt.ext,opt.datapath);
    % user supplied file type
    files = dir(opt.datapath);
    files = {files.name};
    if ~isempty(files)
      for i = 1:length(files)
        [pathstr,name,ext,versn] = fileparts(files{i});
        if strcmp(ext,['.' opt.ext])
          datafile{end+1} = files{i};
        end
      end
    end
  else
    if isfield(opt,'indata')
      % mat datafile
      fprintf('%s: searching for %s in %s\n',mfilename,opt.indata,opt.datapath);
      files = what(opt.datapath);
      if ~isempty(files), files = files.mat; end
      % datafile: check timesurfer data type
      for i=1:length(files)
        if any(strcmp(who('-file',fullfile(opt.datapath,files{i})),opt.indata))
          datafile{end+1} = fullfile(opt.datapath,files{i});
        end
      end
    end
  end
end
opt.datafile = datafile;    

for f = 1:length(opt.datafile)
  if ~exist(opt.datafile{f})
    opt.datafile{f} = fullfile(opt.datapath,opt.datafile{f});
  end
end
opt.datafile = opt.datafile(~strcmp(opt.datafile,'.') & ~strcmp(opt.datafile,'..'));

if isempty(opt.datafile) 
  fprintf('no files found in %s\n',opt.datapath);
  if prompt
    title = sprintf('timesurfer data not found.  Select file.\n');
   	[opt.datafile,opt.datapath] = uigetfile('*.mat',title);
  end
elseif ~isempty(opt.datafile)
  fprintf('%d file(s) found in %s\n',length(opt.datafile),opt.datapath);
end
