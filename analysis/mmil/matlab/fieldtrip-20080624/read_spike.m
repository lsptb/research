function [spike] = read_spike(filename, varargin);

% READ_SPIKE reads spike timestamps and waveforms from various data
% formats.
%
% Use as
%  [spike] = read_spike(filename, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'spikeformat'
%
% The output spike structure contains
%   spike.label     = 1xNchans cell-array, with channel labels
%   spike.waveform  = 1xNchans cell-array, each element contains a matrix (Nsamples X Nspikes)
%   spike.timestamp = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%   spike.unit      = 1xNchans cell-array, each element contains a vector (1 X Nspikes)
%
% See also READ_HEADER, READ_DATA, READ_EVENT

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: read_spike.m,v $
% Revision 1.11  2008/03/25 10:59:02  roboos
% use either NLX_Base_Class_Name or AcqEntName, whichever is available
%
% Revision 1.10  2008/03/04 11:17:07  roboos
% added support for neuralynx_nst (tested) and neuralynx_ntt (untested)
%
% Revision 1.9  2007/03/26 12:32:41  roboos
% changed the API for plexon_plx
%
% Revision 1.8  2007/03/21 13:00:02  roboos
% keep the original data header in the output structure
%
% Revision 1.7  2007/03/19 17:08:57  roboos
% implemented neuralynx_nte
% use timestamp_plexon as low level function instead of replicating the typecasting here
%
% Revision 1.6  2007/03/18 22:02:41  roboos
% also deal with plexon plx spike channels that do not contain any data
%
% Revision 1.5  2007/03/13 14:32:47  roboos
% removed header as optional argument, since read_header does not support spike-only files
% in case of plexon_plx, read the header from the file using low-level importer
% implemented support for plexon_nex, type 0 and 3
%
% Revision 1.4  2007/02/27 09:56:28  roboos
% added some documentation
%
% Revision 1.3  2007/01/09 09:40:38  roboos
% added neuralynx_nse
%
% Revision 1.2  2007/01/04 17:14:12  roboos
% deblank channel labels, renamed data to waveform
%
% Revision 1.1  2007/01/04 12:10:14  roboos
% new implementation, sofar only for plexon_plx
%

% get the options
spikeformat   = keyval('spikeformat',   varargin);

% determine the filetype
if isempty(spikeformat)
  spikeformat = filetype(filename);
end

switch spikeformat
  case {'neuralynx_ncs' 'plexon_ddt'}
    % these files only contain continuous data
    error('file does not contain spike timestamps or waveforms');

  case 'matlab'
    % plain matlab file with a single variable in it
    load(filename, 'spike');

  case 'neuralynx_nse'
    % single channel file, read all records
    nse = read_neuralynx_nse(filename);
    if isfield(nse.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nse.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nse.hdr.AcqEntName};
    end
    spike.timestamp = {nse.TimeStamp};
    spike.waveform  = {nse.dat};
    spike.unit      = {nse.CellNumber};
    spike.hdr       = nse.hdr;

  case 'neuralynx_nst'
    % single channel stereotrode file, read all records
    nst = read_neuralynx_nst(filename, 1, inf);
    if isfield(nst.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nst.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nst.hdr.AcqEntName};
    end
    spike.timestamp = {nst.TimeStamp};
    spike.waveform  = {nst.dat};
    spike.unit      = {nst.CellNumber};
    spike.hdr       = nst.hdr;

  case 'neuralynx_ntt'
    % single channel stereotrode file, read all records
    ntt = read_neuralynx_ntt(filename);
    if isfield(ntt.hdr, 'NLX_Base_Class_Name')
      spike.label   = {ntt.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {ntt.hdr.AcqEntName};
    end
    spike.timestamp = {ntt.TimeStamp};
    spike.waveform  = {ntt.dat};
    spike.unit      = {ntt.CellNumber};
    spike.hdr       = ntt.hdr;

  case 'neuralynx_nte'
    % single channel file, read all records
    nte = read_neuralynx_nte(filename);
    if isfield(nte.hdr, 'NLX_Base_Class_Name')
      spike.label   = {nte.hdr.NLX_Base_Class_Name};
    else
      spike.label   = {nte.hdr.AcqEntName};
    end
    spike.timestamp = {nte.TimeStamp(:)'};
    spike.waveform  = {zeros(0,length(nte.TimeStamp))};  % does not contain waveforms
    spike.unit      = {zeros(0,length(nte.TimeStamp))};  % does not contain units
    spike.hdr       = nte.hdr;

  case 'plexon_nex'
    % a single file can contain multiple channels of different types
    hdr  = read_plexon_nex(filename);
    typ  = [hdr.VarHeader.Type];
    chan = 0;

    spike.label     = {};
    spike.waveform  = {};
    spike.unit      = {};
    spike.timestamp = {};

    for i=1:length(typ)
      if typ(i)==0
        % neurons, only timestamps
        nex = read_plexon_nex(filename, 'channel', i);
        nspike = length(nex.ts);
        chan = chan + 1;
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = zeros(0, nspike);
        spike.unit{chan}      = nan*ones(1,nspike);
        spike.timestamp{chan} = nex.ts;
      elseif typ(i)==3
        % neurons, timestamps and waveforms
        nex = read_plexon_nex(filename, 'channel', i);
        chan = chan + 1;
        nspike = length(nex.ts);
        spike.label{chan}     = deblank(hdr.VarHeader(i).Name);
        spike.waveform{chan}  = nex.dat;
        spike.unit{chan}      = nan*ones(1,nspike);
        spike.timestamp{chan} = nex.ts;
      end
    end
    spike.hdr = hdr;

  case 'plexon_plx'
    % read the header information
    hdr   = read_plexon_plx(filename);
    nchan = length(hdr.ChannelHeader);
    typ   = [hdr.DataBlockHeader.Type];
    unit  = [hdr.DataBlockHeader.Unit];
    chan  = [hdr.DataBlockHeader.Channel];

    for i=1:nchan
      % select the data blocks that contain spike waveforms and that belong to this channel
      sel = (typ==1 & chan==hdr.ChannelHeader(i).Channel);

      if any(sel)
        % get the timestamps that correspond with this spike channel
        tsl = [hdr.DataBlockHeader(sel).TimeStamp];
        tsh = [hdr.DataBlockHeader(sel).UpperByteOf5ByteTimestamp];
        % convert the 16 bit high timestamp into a 32 bit integer
        ts = timestamp_plexon(tsl, tsh);
        spike.timestamp{i} = ts;
        spike.unit{i}      = unit(sel);
      else
        % this spike channel is empty
        spike.timestamp{i} = [];
        spike.unit{i}      = [];
      end
    end
    for i=1:nchan
      spike.label{i}    = deblank(hdr.ChannelHeader(i).Name);
      spike.waveform{i} = read_plexon_plx(filename, 'ChannelIndex', i, 'header', hdr);
    end
    spike.hdr = hdr;

  otherwise
    error('unsupported data format');
end

