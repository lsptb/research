function auto_script(varargin)

in = mmil_args2parms(varargin, { ...
  'options',[],[],...
  'function','',[],...
  'cmdstr','',[],...
  'outfile','',[],...
  'outpath','',[],...
  'fileid','',[],...
  'filename','',[],...
});

opt = getfield(in,'options');
if ~isstruct(opt), return; end
if isnumeric(in.fileid), in.fileid = num2str(in.fileid); end

if (~isfield(opt,'indata')  || isempty(opt.indata))  && (isfield(opt,'itype') && ~isempty(opt.itype))
  opt.indata = opt.itype;
end
if (~isfield(opt,'outdata') || isempty(opt.outdata)) && (isfield(opt,'otype') && ~isempty(opt.otype))
  opt.outdata = opt.otype;
end
  
try
  timesurfer_flag = isfield(opt,'datafile') && ~isempty(opt.datafile) && ...
                  any(opt.indata) && any(opt.outdata); 
catch
  timesurfer_flag = 0;
end

if ~isempty(in.outfile)
  [pathstr name ext versn] = fileparts(in.outfile);
  outfile = in.outfile;
  if ~exist(pathstr,'dir'),
      fprintf('making directory: %s\n',pathstr);
      unix(['mkdir -p ' pathstr]);
  end      
else
  if ~isempty(in.outpath)
    outfile = in.outpath;
    if ~exist(outfile,'dir'),
        fprintf('making directory: %s\n',outfile);
        unix(['mkdir -p ' outfile]);
    end            
  elseif isfield(opt,'rootoutdir')
    outfile = fullfile(opt.rootoutdir,'scripts');
    if ~exist(outfile,'dir'),
        fprintf('making scripts directory: %s\n',outfile);
        unix(['mkdir -p ' outfile]);
    end        
  else
    outfile = pwd;
  end
	if isfield(opt,'prefix')
		prefix = opt.prefix;
	else
		prefix = 'run';
	end
  if ~isempty(in.function)
    outfile = sprintf('%s/%s_%s%s.m',outfile,prefix,in.function,in.fileid);
  else
    outfile = sprintf('%s/%s_auto_options.m',outfile,prefix);
  end    
end
try opt = rmfield(opt,'previous'); end
fname = fieldnames(opt);

fprintf('%s: writing auto script: %s\n',mfilename,outfile);
fid = fopen(outfile,'wt'); 
fprintf(fid,'tic\n');

if isfield(opt,'mmilclusterheader')
  fprintf(fid,'%s\n',opt.mmilclusterheader);
end

% if timesurfer_flag
if any(opt.indata) && ~isempty(opt.datafile)
  if ~iscell(opt.datafile), opt.datafile = {opt.datafile}; end
  for i = 1:length(opt.datafile)
    fprintf(fid,'datafile{%d} = ''%s'';\n',i,opt.datafile{i});
  end
  fprintf(fid,'for i=1:length(datafile)\n');
  fprintf(fid,'\t%s(i) = getfield(load(datafile{i},''%s''),''%s'');\n',opt.indata,opt.indata,opt.indata);
  fprintf(fid,'end\n');
  fprintf(fid,'if i>1, %s = ts_combine_data(%s); end\n\n',opt.indata,opt.indata);
end

for i = 1:length(fname)
    str = fname{i};
    val = opt.(str);
    valstr = '';
    if isnumeric(val)
        if length(val)==1
            valstr = num2str(val);
        else
            valstr = numarray2str(val);
        end
    elseif iscell(val)
        valstr = '';
        for j = 1:length(val)
            if isstr(val{j})
                if j==1
                    valstr=sprintf('''%s''',val{j}); 
                else
                    valstr = sprintf('%s,''%s''',valstr,val{j});
                end
            elseif isnumeric(val{j})
                if j==1, valstr=numarray2str(val{j}); else valstr=[valstr ',' numarray2str(val{j})]; end
            end
        end
        valstr = sprintf('{%s}',valstr);
    elseif isstruct(val)
        valstr = str;
    elseif isstr(val)
        valstr = sprintf('''%s''',val);
    end
    fprintf(fid,'parms.%s = %s;\n',str,valstr);
end

if timesurfer_flag
  fprintf(fid,'\nargs = mmil_parms2args(parms);\n');
  fprintf(fid,'%s = %s(%s,args{:});\n',opt.outdata,opt.function,opt.indata);
%   fprintf(fid,'try opt.previous = %s.opt; end\n',opt.indata);
%   fprintf(fid,'%s.opt = opt;\n',opt.outdata);
%   fprintf(fid,'save(''%s'',''%s'');\n',in.filename,opt.outdata);
elseif ~isempty(in.function)
	fprintf(fid,'\nargs = mmil_parms2args(parms);\n');
	if ~isempty(in.cmdstr)
		fprintf(fid,'%s\n',in.cmdstr);
	else
  	fprintf(fid,'%s(args{:})\n',in.function);
	end
end
if any(opt.outdata)
  if any(opt.indata), fprintf(fid,'try parms.previous = %s.parms; end\n',opt.indata); end
  fprintf(fid,'%s.parms = parms;\n',opt.outdata);  
  fprintf(fid,'save(''%s'',''%s'');\n',in.filename,opt.outdata);
end
fprintf(fid,'\ntoc\nexit\n');  % may need for cluster computing
fclose(fid);


function str = numarray2str(array)
str = '';
for i = 1:length(array)
    str = sprintf('%s %g',str,array(i));
end
str = sprintf('[%s]',str);

