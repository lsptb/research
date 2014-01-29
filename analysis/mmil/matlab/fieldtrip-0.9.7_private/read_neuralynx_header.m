function [hdr] = read_neuralynx_hdr(dirname);

% READ_NEURALYNX_HDR reads the header information of multiple continuous
% datafiles in a directroy
%
% A Neuralynx dataset consists of separate files, one for each channel. All
% Neuralynx datafiles starts with a 16k header (in ascii format), followed
% by an arbitrary number of data records. The format of the data records
% depend on the type of data contained in the channel (e.g. continuous or
% spike data).
%
% Use as
%   hdr = read_neuralynx_hdr(datadir)
%
% See also READ_NEURALYNX_DATA

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: read_neuralynx_header.m,v $
% Revision 1.6  2005/12/02 09:01:11  roboos
% only read timestamps if file is not empty
%
% Revision 1.5  2005/09/09 08:38:36  roboos
% added the last timestamp and the number of records to the header
%
% Revision 1.4  2005/09/08 16:53:20  roboos
% fixed small error in help message
%
% Revision 1.3  2005/08/03 10:38:40  roboos
% fixed sampling frequency for our Digital Lynx amplifier
%
% Revision 1.2  2005/05/19 07:06:25  roboos
% added FirstTimeStamp to the header for each channel
%
% Revision 1.1  2005/05/09 11:39:38  roboos
% new implementation based on Neuralynx documentation on continuous data files
%

file = dir(dirname);
file = file(~cell2mat({file.isdir}));
for i=1:length(file)
  filename{i} = fullfile(dirname, file(i).name);
end
clear file

for i=1:length(filename)
  ncs(i) = filetype(filename{i}, 'neuralynx_ncs');
end
filename = filename(find(ncs));

orig = {};
for i=1:length(filename)
  orig{i} = getheader(filename{i});
  NumRecords(i)     = numrecords(filename{i});
  if NumRecords(i)>0
    FirstTimeStamp(i) = timestamp(filename{i}, 1);
    LastTimeStamp(i)  = timestamp(filename{i}, inf);
  else
    warning(sprintf('the file ''%s'' does not contain any data', filename{i}));
    FirstTimeStamp(i) = 0;
    LastTimeStamp(i)  = 0;
  end
end

for i=1:length(filename)
  hdr.label{i}              = orig{i}.NLX_Base_Class_Name;
  hdr.SamplingFrequency(i)  = orig{i}.SamplingFrequency;
  hdr.ADBitVolts(i)         = orig{i}.ADBitVolts;
  hdr.FirstTimeStamp        = FirstTimeStamp;
  hdr.LastTimeStamp         = LastTimeStamp;
  hdr.NumRecords            = NumRecords;
end

if any(hdr.SamplingFrequency~=hdr.SamplingFrequency(1))
  warning('not all channels have the same sampling rate');
  hdr.SamplingFrequency = median(hdr.SamplingFrequency);  % FIXME, this is silly
else
  hdr.SamplingFrequency = hdr.SamplingFrequency(1);
end

warning('fixing sampling frequency to 32556 Hz');         % FIXME, this is even more silly
hdr.SamplingFrequency = 32556;

% remember the detailled header information for each channel
hdr.orig = orig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading the header of a single channel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = getheader(filename);
fid     = fopen(filename, 'rb', 'ieee-le');
buf     = fread(fid, 16*1024, 'char');
buf     = buf(:)';
nl      = find(buf==10);    % determine the new-lines
cr      = find(buf==13);    % determine the carriage-returns
begchar = [1 nl(1:(end-1))];
endchar = nl - 1;
num     = length(nl);
hdr     = [];

for i=1:num
  line = fliplr(deblank(fliplr(deblank(char(buf(begchar(i):endchar(i)))))));
  if line(1)=='#'
    % line contains a comment
    continue
  else
    % strip the '-' sign
    if line(1)=='-'
      line = line(2:end);
    end
    % replace tabs with spaces
    line(find(line==9)) = ' ';
    % cut into pieces
    item = strread(line, '%s');
    if length(item)==2
      key = item{1};
      val = item{2};
      if any(val(1)=='-01234567989')
        % try to convert to number
        val = str2num(val);
        if isempty(val)
          % revert to the original text
          val = item{2};
        end
      end
      % assign the value to the header structure
      hdr = setfield(hdr, key, val);
    end
  end
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading a single timestamp of a single channel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t] = timestamp(filename, num)
fid = fopen(filename, 'rb', 'ieee-le');
headersize = 16*1024;
recordsize = 8+3*4+512*2;
if ~isinf(num)
  % read the timestamp of the indicated record
  fseek(fid, headersize + (num-1)*recordsize, 'bof');
  t = fread(fid, 1, 'int64');
else
  % read the timestamp of the last record
  fseek(fid, -recordsize, 'eof');
  t = fread(fid, 1, 'int64');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for determining the number of records in a single channel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t] = numrecords(filename);
fid = fopen(filename, 'rb', 'ieee-le');
headersize = 16*1024;
recordsize = 8+3*4+512*2;
fseek(fid, 0, 'eof');
t = (ftell(fid) - headersize)/recordsize;
fclose(fid);

