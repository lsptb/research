function odat = process_protocol(protocol,varargin)
% get datafiles from
% if data_source is a number (index to input step)
%   protocol.steps(data_source).parms.filename
%   protocol.steps(data_source).parms.hdr
%     load hdr => hdr.parms.filename
% elseif data_source is a string containing 'hdr' (hdr filename)
%   load hdr => hdr.parms.filename
% elseif data_source is a string
%   data_source
% elseif data_source is empty
%   protocol.steps(this).parms.datafile
%   protocol.steps(this).parms.hdr

inparms = mmil_args2parms(varargin,...
						{'log',0,{0,1},...              
						 'run',1,{0,1},...
						},false);

% get protocol structure
if ischar(protocol) && strfind(protocol,'.csv')
  protocol = read_protocol(protocol);
elseif ischar(protocol) && strfind(protocol,'.mat')
  load(protocol);
elseif ~isstruct(protocol)
  error('%s: protocol must be passed as a structure or csv or mat filename',mfilename);
end

% get steps to process
s    = protocol.run{end};

% process steps
for n = 1:length(s)
  % index to this stepid
  k = s(n);
  % info for this step
  step  = protocol.steps(k);  % protocol info
  fparm = step.parms;         % funparms for run_timesurfer
  parms = fparm.parms;        % parms for TS functions
  fprintf('%s: processing function %g of %g in protocol: %s\n',mfilename,n,length(s),step.function);
  if ischar(step.parms.itype) && exist(step.parms.itype,'var') && step.data_source==s(n-1)
    eval(fparm.funcall);
    continue;
  end
    % do nothing
  % set datafile or header filename
  if isfield(parms,'datafile') && ischar(parms.datafile)
    % do nothing
  elseif isfield(parms,'hdr') && ischar(parms.hdr)
    % do nothing
  elseif isnumeric(step.data_source)
    % get header or datafile from the 1st data_source index
    datidx = [protocol.steps.stepid] == step.data_source(1);
    if isfield(protocol.steps(datidx),'hdr')
      parms.hdr = protocol.steps(datidx).hdr;
    elseif issubfield(protocol.steps(datidx).parms,'parms.filename')
      parms.datafile = protocol.steps(datidx).parms.parms.filename;
    end
  elseif ischar(step.data_source) && strfind(step.data_source,'hdr')
    % set as header
    parms.hdr = step.data_source;
  elseif ischar(step.data_source)
    % this is a datafile
    parms.datafile = step.data_source;
  end
  % CALL FUNCTION
  fparm.specflags.save_flag = step.save_flag;
  fparm.parms = parms;
  clear parms;
  parms.function = fparm;
  if step.save_flag == 3
    eval(sprintf('%s = run_timesurfer(parms);',fparm.otype));
  else
    [hdrfiles datfiles] = run_timesurfer(fparm);
    protocol.steps(k).hdr = hdrfiles;
    protocol.steps(k).filename = datfiles;
  end
  clear parms fparm step
end

    
    