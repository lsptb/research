function protocol = read_protocol(csvfile)

[data, result]  = mmil_readtext(csvfile, ',','','','empty2zero');

rows = find(result.numberMask(:,1));
cnt  = 1;
for r = 1:length(rows)
  this = data(rows(r),:);
  pres = length(this) >= 1:8;
  if pres(1), stepid  = this{1}; else continue; end
  if pres(2), funname = this{2}; else funname = sprintf('fun%g',num2str(cnt)); end
  if pres(3), parmsrc = this{3}; else parmsrc = ''; end
  if pres(4), datasrc = this{4}; else datasrc = ''; end
  if pres(5), procflg = this{5}; else procflg = 0;  end
  if pres(6), saveflg = this{6}; else saveflg = 0;  end
  if pres(7), outtype = this{7}; else outtype = ''; end
  if pres(8), pause   = this{8}; else pause   = 0;  end
  protocol.steps(cnt).stepid        = stepid;   % unique id for this step
  protocol.steps(cnt).function      = funname;  % name of function executed this step
  protocol.steps(cnt).parm_source   = parmsrc;  % source of function-specific parameters
  if isnumeric(datasrc)
    % datasrc is an index to the input step
    protocol.steps(cnt).data_source = datasrc;
  elseif ~isempty(str2num(datasrc))
    % datasrc is a vector of indices to input steps
    protocol.steps(cnt).data_source = str2num(datasrc);
  elseif ischar(datasrc)
    % datasrc is a filename
    protocol.steps(cnt).data_source = datasrc;
  end
  protocol.steps(cnt).proc_flag     = procflg;  % whether to run this function
  protocol.steps(cnt).pause         = pause;    % whether it needs to wait on a prior function
  protocol.steps(cnt).save_flag     = saveflg;  % whether to save the output
  protocol.steps(cnt).otype         = outtype;  % optional specification of the output to save
  if ~isfield(protocol.steps(cnt),'parms')      % set the function-specific parameters
    protocol.steps(cnt).parms       = [];
  end
  % load function parameters from csv or matfile if procflag ~= 0
  if procflg && exist(parmsrc,'file') && isempty(protocol.steps(cnt).parms)
    if strfind(parmsrc,'.csv')
      parms = load_parms(parmsrc);
      fidx  = ismember({parms.function.name},funname);
      protocol.steps(cnt).parms = parms.function(fidx);
    elseif strfind(parmsrc,'.mat')
      tmp = load(parmsrc);
      fld = fieldnames(tmp);
      if isfield(tmp.(fld{1}),'parms')
        protocol.steps(cnt).parms = tmp.(fld{1}).parms;
      end
      clear tmp
    end
  end 
  cnt = cnt + 1;
end
if ~isfield(protocol,'run')
  protocol.run{1} = find([protocol.steps.proc_flag]);
else
  protocol.run{end+1} = find([protocol.steps.proc_flag]);
end
[path name] = fileparts(csvfile);
protocol.parms.filename{1} = [name '.mat'];
if nargout == 0
  save(protocol.parms.filename{1},'protocol')
end
