function [dat] = read_neuralynx_data(dirname, hdr, begsample, endsample, chanindx);

% READ_NEURALYNX_DATA reads data from multiple continuous datafiles
% that are stored together in a directory.
% 
% Use as 
%   [dat] = read_neuralynx_data(datadir, hdr, begsample, endsample, chanindx)
% where
%   datadir   = directory name
%   hdr       = structure with header information
%   begsample = first sample to read
%   endsample = last sample to read
%   chanindx  = list of channel numbers to read
%
% Samples are counted starting at 1, which corresponds to the first 
% sample that is present for all channels.
%
% Channels are counted according to their order in the common header
% structure.
%
% The output signal is expressed in microvolts (uV).
% 
% See also READ_NEURALYNX_HEADER

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: read_neuralynx_data.m,v $
% Revision 1.12  2006/02/15 09:25:03  roboos
% once more: complete redesign of the reading
% current version gives warning in case of discontinuous data (should be an error)
% current version gives warning in case of different FirstTimeStamp over channels (should be an error)
%
% Revision 1.11  2005/12/02 09:00:27  roboos
% made the copying-pasting of the data segments more efficient
%
% Revision 1.10  2005/11/23 15:49:49  roboos
% fixed bug in selection of desired samples if begrecord~=1
%
% Revision 1.9  2005/09/09 12:27:59  roboos
% change in help: output unit (uV) is explicitely mentioned
%
% Revision 1.8  2005/09/09 08:40:32  roboos
% fixed a bug in the selection of samples from the time-shifted trials
%
% Revision 1.7  2005/08/05 13:40:02  roboos
% removed 1000x scaling, since that was due to the input mouse
%
% Revision 1.6  2005/08/04 07:45:26  roboos
% fixed a bug in the channel numbering
%
% Revision 1.5  2005/06/23 15:30:42  roboos
% fixed bug in begin and endsamples of the segment of interest after reading in the buffer
%
% Revision 1.4  2005/06/01 08:04:46  roboos
% added comment in code
%
% Revision 1.3  2005/05/19 07:05:48  roboos
% moved reading of timestamps to separate subfunction
%
% Revision 1.2  2005/05/12 07:23:41  roboos
% added check for trying to read beyond the end of file
%
% Revision 1.1  2005/05/09 11:39:38  roboos
% new implementation based on Neuralynx documentation on continuous data files
%

% determine all the files that are contained in this dataset
file = dir(dirname);
file = file(~cell2mat({file.isdir}));
for i=1:length(file)
  filename{i} = fullfile(dirname, file(i).name);
end
clear file

% only select the continuous channels
for i=1:length(filename)
  ncs(i) = filetype(filename{i}, 'neuralynx_ncs');
end
filename = filename(find(ncs));

% the default is to read all continuous channels
if nargin<5
  chanindx = 1:length(filename);
end
% it should be a row vector
chanindx = chanindx(:)';

% the channel ordering should correspond with the labels in the header structure
for i=1:length(chanindx)
  thischan = chanindx(i);
  [p, f, x]  = fileparts(filename{thischan});
  if ~strcmp(hdr.label{thischan}, f)
    error('channel names in the header and the data directory do not match');
  end
end

% The file starts with a 16*1024 bytes header in ascii, followed by a
% number of records (c.f. trials).
%
% The format of a continuous sampled record is
%   int64 TimeStamp
%   int32 ChanNumber
%   int32 SampFreq
%   int32 NumValidSamp
%   int16 Samp[0] ... int16 Samp[511]
% Note that if NumValidSamp < 512, Samp[n], where n >= NumValidSamp, will
% contain random data.

% work with the original header, not with the FieldTrip header
dum = hdr;
hdr = hdr.orig;

% one timestamp is 1us, and one sample is 1/hdr.Fs, however there is a slight inaccuracy: the 32556Hz seems to be 32552Hz
% TimeStampPerSample = 1e6/hdr.Fs;
% determine the exact number of timestamps per sample
TimeStampPerRecord = (hdr.LastTimeStamp - hdr.FirstTimeStamp)./(hdr.NumRecords-1);
TimeStampPerSample = TimeStampPerRecord/512;

if any(TimeStampPerSample~=TimeStampPerSample(1))
  % error('the data is not continuous in all channels');
  warning('the data is not continuous in all channels');
end
TimeStampPerSample = TimeStampPerSample(1);

if any(hdr.FirstTimeStamp~=hdr.FirstTimeStamp(1))
  % error('the data does not start at the same sample in all channels');
  warning('the data does not start at the same sample in all channels');
end

nsamples = endsample-begsample+1;
nchans   = length(chanindx);
dat      = zeros(nchans, nsamples);

% include a scaling factor into ADBitVolts to convert to uV
hdr.ADBitVolts = 1e6 * hdr.ADBitVolts;
% note that the so-called "mouse" that connects to the HS-27 headstage has an attenuation of 1000x
% an external calibration signal therefore appears to be 1000x too small

% determine the records that contain the sample numbers of the requested segment
begrecord  = floor(begsample/512) + 1;
endrecord  = ceil (endsample/512);
nrecord    = endrecord-begrecord+1;
buf        = zeros(1,nrecord*512);

for i=1:length(chanindx)
  thischan   = chanindx(i);
  fid = fopen(filename{thischan}, 'rb', 'ieee-le');
  fseek(fid, 16*1024, 'bof');                   % skip the header
  fseek(fid, (begrecord-1)*(512*2+20), 'cof');  % skip to the first record
  for t=begrecord:endrecord
    % read a single continuous data record
    TimeStamp    = fread(fid,   1, 'int64');
    ChanNumber   = fread(fid,   1, 'int32');
    SampFreq     = fread(fid,   1, 'int32');
    NumValidSamp = fread(fid,   1, 'int32');
    Samp         = fread(fid, 512, 'int16');
    if length(Samp)<512
      error('cannot read beyond end of file');
    end
    % mark the invalid samples
    Samp((NumValidSamp+1):end) = nan;
    % add this record to the output data
    buf(((t-begrecord)*512+1):(t-begrecord+1)*512) = Samp;
  end
  fclose(fid);
  % calibrate the segment of data that is contained in this record 
  % ADBitVolts also contains a scaling factor to convert to uV
  buf = hdr.ADBitVolts(thischan)*buf;
  % select the desired samples out of the buffer
  begsel = mod(begsample, 512); 
  endsel = begsel + nsamples - 1;
  dat(i,:) = buf(begsel:endsel);
end

