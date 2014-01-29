function [trialData] = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx)

% READ_SBIN_DATA reads the data from an EGI segmented simple binary format file
%
% Use as
%   [trialData] = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx)
% with
%   filename       name of the input file
%   hdr            header structure, see READ_HEADER
%   begtrial       first trial to read, mutually exclusive with begsample+endsample
%   endtrial       last trial to read,  mutually exclusive with begsample+endsample
%   chanindx       list with channel indices to read
%
% This function returns a 3-D matrix of size Nchans*Nsamples*Ntrials.
%_______________________________________________________________________
%
%
% Modified from EGI's readEGLY.m with permission 2008-03-31 Joseph Dien
%

fh=fopen([filename],'r');
if fh==-1
  error('wrong filename')
end

version		= hdr.orig.header_array(1);
precision = bitand(version,6);
Nevents=hdr.orig.header_array(17);

switch precision
    case 2
        trialLength=2*hdr.nSamples*(hdr.nChans+Nevents)+6;
        dataType='int16';
    case 4
        trialLength=4*hdr.nSamples*(hdr.nChans+Nevents)+6;
        dataType='single';
    case 6
        trialLength=8*hdr.nSamples*(hdr.nChans+Nevents)+6;
        dataType='double';
end

fseek(fh, 40+length(hdr.orig.CatLengths)+sum(hdr.orig.CatLengths)+Nevents*4, 'cof'); %skip over header
fseek(fh, (begtrial-1)*trialLength, 'cof'); %skip over initial segments

trialData=zeros(hdr.nChans,hdr.nSamples,endtrial-begtrial+1);

for segment=1:(endtrial-begtrial+1)
    fseek(fh, 6, 'cof'); %skip over segment info
    temp = fread(fh, [(hdr.nChans+Nevents), hdr.nSamples],dataType);
    trialData(:,:,segment) = temp(1:hdr.nChans,:);
end
trialData=trialData(chanindx, :,:);

fclose(fh);