% pop_loadbv() - load Brain Vision Data Exchange format dataset and
%                return EEGLAB EEG structure
%
% Usage:
%   >> [EEG, com] = pop_loadbv; % pop-up window mode
%   >> [EEG, com] = pop_loadbv(path, hdrfile);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, srange);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, [], chans);
%   >> [EEG, com] = pop_loadbv(path, hdrfile, srange, chans);
%
% Optional inputs:
%   path      - path to files
%   hdrfile   - name of Brain Vision vhdr-file (incl. extension)
%   srange    - scalar first sample to read (up to end of file) or
%               vector first and last sample to read (e.g., [7 42];
%               default: all)
%   chans     - vector channels channels to read (e.g., [1:2 4];
%               default: all)
%
% Outputs:
%   EEG       - EEGLAB EEG structure
%   com       - history string
%
% Author: Andreas Widmann & Arnaud Delorme, 2004-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: pop_loadbv.m,v $
% Revision 1.11  2005/12/19 16:12:24  andreaswidmann
% Readded EEG.event.urevent field
%
% Revision 1.10  2005/12/16 13:11:29  andreaswidmann
% Removed EEG.chanlocs.type field
%
% Revision 1.9  2005/12/13 11:34:42  andreaswidmann
% Unmatched end with MATLAB < R14
%
% Revision 1.8  2005/11/29 16:19:34  andreaswidmann
% Typo in memory preallocation, removed EEG.chanlocs.datachan field
%
% Revision 1.7  2005/11/24 13:00:42  andreaswidmann
% Marker duration
%
% Revision 1.6  2005/11/22 18:59:57  arnodelorme
% Updating pop_loadbv & pop_writebva
%
% Revision 1.5  2005/11/17 17:16:56  andreaswidmann
% MATLAB < R14 compatibility issues
%
% Revision 1.4  2005/11/08 14:30:10  andreaswidmann
% Log
%
% Revision 1.3  2005/11/08 09:22:34  andreaswidmann
% Double output precision for MATLAB < R14
%
% Revision 1.2  2005/11/07 10:31:01  andreaswidmann
% Parameter plausibility checking
%
% Revision 1.1.1.1  2005/11/03 22:57:36  arnodelorme
% initial import into CVS
%
% Revision 1.5  2005/11/05 19:30:00  Stefan Debener
% Readconfig backward compatibility
%
% Revision 1.4  2005/11/02 15:12:00 andreaswidmann
% Coordinates, added srange parameter, merge revisions 1.1.1.1 and 1.3,
% memory and speed optimization of data file reading strategies,
% SegmentationType, GPL, config file parsing rewritten
%
% Revision 1.3  2005/09/14 08:21:00  arnodelorme
% DataOrientation, BinaryFormat, homogenous boundary intervals, removed
% srange parameter
%
% Revision 1.2  2005/05/11 18:42:00  andreaswidmann
% Parsing of channel info (fileformat change in Vision Recorder 1.03)
%
% Revision 1.1.1.1 2005/08/05 18:45:00  msiegel
% Added chans parameter
%
% Revision 1.1  2004/12/21 15:19:00  andreaswidmann
% Added srange parameter
%
% Revision 1.0.1.1  2004/11/17 02:53:00  arnodelorme
% Formating for EEGLAB
%
% Revision 1.0  2004/11/16 15:11:00  andreaswidmann
% Initial revision
%

function [EEG, com] = pop_loadbv(path, hdrfile, srange, chans)

    com = '';
    EEG = [];

    if nargin < 2
        [hdrfile path] = uigetfile2('*.vhdr', 'Select Brain Vision vhdr-file - pop_loadbv()');
        if hdrfile(1) == 0, return; end

        drawnow;
        uigeom = {[1 0.5] [1 0.5]};
        uilist = {{ 'style' 'text' 'string' 'Interval (samples; e.g., [7 42]; default: all):'} ...
                  { 'style' 'edit' 'string' ''} ...
                  { 'style' 'text' 'string' 'Channels (e.g., [1:2 4]; default: all):'} ...
                  { 'style' 'edit' 'string' ''}};
        result = inputgui(uigeom, uilist, 'pophelp(''pop_loadbv'')', 'Load a Brain Vision Data Exchange format dataset');
        if length( result ) == 0, return; end
        if ~isempty(result{1}),
            srange = str2num(result{1});
        end
        if ~isempty(result{2}),
            chans = str2num(result{2});
        end
    end

    % Header file
    disp('pop_loadbv(): reading header file');
    hdr = readconfig(path, hdrfile);

    % Common Infos
    try, EEG = eeg_emptyset;
    catch, end;
    EEG.comments = ['Original file: ' hdr.commoninfos.datafile];
    if ~strcmpi(hdr.commoninfos.dataformat, 'binary')
        error('Ascii data file import not (yet?) implemented.');
    end
    hdr.commoninfos.numberofchannels = str2num(hdr.commoninfos.numberofchannels);
    EEG.srate = 1000000 / str2num(hdr.commoninfos.samplinginterval);

    % Binary Infos
    switch lower(hdr.binaryinfos.binaryformat)
        case 'int_16',        binformat = 'int16'; bps = 2;
        case 'uint_16',       binformat = 'uint16'; bps = 2;
        case 'ieee_float_32', binformat = 'float32'; bps = 4;
        otherwise error('Unsupported binary format');
    end

    % Channel Infos
    if ~exist('chans', 'var') || isempty(chans)
        chans = 1:hdr.commoninfos.numberofchannels;
        EEG.nbchan = hdr.commoninfos.numberofchannels;
    else
        EEG.nbchan = length(chans);
    end
    if any(chans < 1) || any(chans > hdr.commoninfos.numberofchannels)
        error('chans out of available channel range');
    end
    if isfield(hdr, 'channelinfos')
        for chan = 1:length(chans)
            [EEG.chanlocs(chan).labels, chanlocs(chan).scale] = strread(hdr.channelinfos{chans(chan)}, '%s%*s%f', 1, 'delimiter', ',');
            EEG.chanlocs(chan).labels = char(EEG.chanlocs(chan).labels);
%             EEG.chanlocs(chan).datachan = chans(chan);
        end
        if isempty([chanlocs.scale])
            chanlocs = rmfield(chanlocs, 'scale');
        end
    end;
%     [EEG.chanlocs.type] = deal([]);

    % Coordinates
    if isfield(hdr, 'coordinates')
        for chan = 1:length(chans)
            [EEG.chanlocs(chan).sph_radius, theta, phi] = strread(hdr.coordinates{chans(chan)}, '%f%f%f', 'delimiter', ',');
            if EEG.chanlocs(chan).sph_radius == 0 && theta == 0 && phi == 0
                EEG.chanlocs(chan).sph_radius = [];
                EEG.chanlocs(chan).sph_theta = [];
                EEG.chanlocs(chan).sph_phi = [];
            else
                EEG.chanlocs(chan).sph_theta = phi - 90 * sign(theta);
                EEG.chanlocs(chan).sph_phi = -abs(theta) + 90;
            end
        end
        try
            [EEG.chanlocs, EEG.chaninfo] = pop_chanedit(EEG.chanlocs, 'convert', 'sph2topo');
            [EEG.chanlocs, EEG.chaninfo] = pop_chanedit(EEG.chanlocs, 'convert', 'sph2cart');
        catch, end
    end

    % Open data file
    disp('pop_loadbv(): reading EEG data');
    [IN, message] = fopen(fullfile(path, hdr.commoninfos.datafile));
    if IN == -1
        [IN, message] = fopen(fullfile(path, lower(hdr.commoninfos.datafile)));
        if IN == -1
            error(message);
        end;
    end
    if isfield(hdr.commoninfos, 'datapoints')
        hdr.commoninfos.datapoints = str2num(hdr.commoninfos.datapoints);
    else
        fseek(IN, 0, 'eof');
        hdr.commoninfos.datapoints = ftell(IN) / (hdr.commoninfos.numberofchannels * bps);
        fseek(IN, 0, 'bof');
    end

    % Sample range
    if ~exist('srange', 'var') || isempty(srange)
        srange = 1;
        EEG.pnts = hdr.commoninfos.datapoints;
    elseif length(srange) == 1
        EEG.pnts = hdr.commoninfos.datapoints - srange(1) + 1;
    else
        EEG.pnts = srange(2) - srange(1) + 1;
    end
    if any(srange < 1) || any(srange > hdr.commoninfos.datapoints)
        error('srange out of available data range');
    end

    % Read data
    switch lower(hdr.commoninfos.dataorientation)
        case 'multiplexed'
            if EEG.nbchan == hdr.commoninfos.numberofchannels % Read all channels
                fseek(IN, (srange(1) - 1) * EEG.nbchan * bps, 'bof');
                EEG.data = fread(IN, [EEG.nbchan, EEG.pnts], [binformat '=>float32']);
            else % Read channel subset
                EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]); % Preallocate memory
                for chan = 1:length(chans)
                    fseek(IN, (srange(1) - 1) * hdr.commoninfos.numberofchannels * bps + (chans(chan) - 1) * bps, 'bof');
                    EEG.data(chan, :) = fread(IN, [1, EEG.pnts], [binformat '=>float32'], (hdr.commoninfos.numberofchannels - 1) * bps);
                end
            end
        case 'vectorized'
            if isequal(EEG.pnts, hdr.commoninfos.datapoints) && EEG.nbchan == hdr.commoninfos.numberofchannels % Read entire file
                EEG.data = fread(IN, [EEG.pnts, EEG.nbchan], [binformat '=>float32']).';
            else % Read fraction of file
                EEG.data = repmat(single(0), [EEG.nbchan, EEG.pnts]); % Preallocate memory
                for chan = 1:length(chans)
                    fseek(IN, ((chans(chan) - 1) * hdr.commoninfos.datapoints + srange(1) - 1) * bps, 'bof');
                    EEG.data(chan, :) = fread(IN, [1, EEG.pnts], [binformat '=>float32']);
                end
            end
        otherwise
            error('Unsupported data orientation')
    end
    fclose(IN);
    EEG.trials = 1;
    EEG.xmin   = 0;
    EEG.xmax   = (EEG.pnts - 1) / EEG.srate;

    % Convert to EEG.data to double for MATLAB < R14
    if str2num(version('-release')) < 14
        EEG.data = double(EEG.data);
    end

    % Scale data
    if exist('chanlocs', 'var') && isfield(chanlocs, 'scale')
        disp('pop_loadbv(): scaling EEG data');
        for chan = 1:EEG.nbchan
            EEG.data(chan, :) = EEG.data(chan, :) * chanlocs(chan).scale;
        end
    end

    % Marker file
    disp('pop_loadbv(): reading marker file');
    mrk = readconfig(path, hdr.commoninfos.markerfile);

    if hdr.commoninfos.datafile ~= mrk.commoninfos.datafile
        disp('pop_loadbv() warning: data files in header and marker files inconsistent.');
    end

    % Marker infos
    if isfield(mrk, 'markerinfos')
        for index = 1:length(mrk.markerinfos)
            [mrktype, mrkdesc, EEG.event(index).latency, mrksizepnts, mrkchan, mrktime] = ...
                strread(mrk.markerinfos{index}, '%s%s%f%d%d%d', 'delimiter', ',');
            if strcmpi(mrktype, 'New Segment') || strcmpi(mrktype, 'DC Correction')
                EEG.event(index).type = 'boundary';
                EEG.event(index).code = char(mrktype);
                EEG.event(index).duration = mrksizepnts;
            else
                EEG.event(index).type = char(mrkdesc);
                EEG.event(index).code = char(mrktype);
                EEG.event(index).duration = mrksizepnts;
            end
            EEG.event(index).latency = EEG.event(index).latency - srange(1) + 1;
            EEG.event(index).urevent = index;
        end

        EEG.event = EEG.event(find([EEG.event.latency] >= 1 & [EEG.event.latency] <= EEG.pnts));
        EEG.urevent = rmfield(EEG.event, 'urevent');

        % find if boundaries at homogenous intervals
        % ------------------------------------------
        boundaries = strmatch('boundary', {EEG.event.type});
        boundlats = unique([EEG.event(boundaries).latency]);
        if (isfield(hdr.commoninfos, 'segmentationtype') && (strcmpi(hdr.commoninfos.segmentationtype, 'markerbased') || strcmpi(hdr.commoninfos.segmentationtype, 'fixtime'))) && length(boundaries) > 1 && length(unique(diff([boundlats EEG.pnts + 1]))) == 1
            EEG.trials = length(boundlats);
            EEG.pnts   = EEG.pnts / EEG.trials;
            EEG.event(boundaries) = [];

            % adding epoch field
            % ------------------
            for index = 1:length(EEG.event)
                EEG.event(index).epoch = ceil(EEG.event(index).latency / EEG.pnts);
            end

            % finding minimum time
            % --------------------
            tles = strmatch('time 0', lower({EEG.event.code}))';
            if ~isempty(tles)
                [EEG.event(tles).type] = deal('TLE');
                EEG.xmin = -(EEG.event(tles(1)).latency - 1) / EEG.srate;
            end
        end
    end

    EEG.ref = 'common';

    try, EEG = eeg_checkset(EEG);
    catch, end

    if nargout == 2
        com = sprintf('EEG = pop_loadbv(''%s'', ''%s'', %s, %s);', path, hdrfile, mat2str(srange), mat2str(chans));
    end

% Read config file
function config = readconfig(path, file)

    if nargin < 2
        error('Not enough input arguments');
    end

    % Open and read file
    [IN, message] = fopen(fullfile(path,file));
    if IN == -1
        [IN, message] = fopen(fullfile(path, lower(file)));
        if IN == -1
            error(message);
        end;
    end
    raw={};
    while ~feof(IN)
        raw = [raw; {fgetl(IN)}];
    end
    fclose(IN);

    % Remove comments and empty lines
    raw(strmatch(';', raw)) = [];
    raw(find(cellfun('isempty', raw) == true)) = [];

    % Find sections
    sections = [strmatch('[', raw)' length(raw) + 1];
    for section = 1:length(sections) - 1

        % Convert section name
        fieldname = lower(char(strread(raw{sections(section)}, '[%s', 'delimiter', ']')));
        fieldname(find(isspace(fieldname) == true)) = [];

        % Fill structure with parameter value pairs
        switch fieldname
            case {'commoninfos' 'binaryinfos'}
                for line = [sections(section) + 1:sections(section + 1) - 1]
                    [parameter, value] = strread(raw{line}, '%s%s', 'delimiter', '=');
                    config.(fieldname).(lower(char(parameter))) = char(value);
                end
            case {'channelinfos' 'coordinates' 'markerinfos'}
                for line = [sections(section) + 1:sections(section + 1) - 1]
                    [parameter, value] = strread(raw{line}, '%s%s', 'delimiter', '=');
                    config.(fieldname)(str2num(parameter{1}(3:end))) = value;
                end
            case 'comment'
                config.(fieldname) = raw(sections(section) + 1:sections(section + 1) - 1);
        end
    end
