function output = importAxoX(fn)
%Import Axograph X files
% DOES NOT WORK ON OLDER AXOGRAPH FILES!

if nargin < 1
    [fn, pn] = uigetfile('*.*', 'Pick an Axograph file', 'MultiSelect','on');
else
    pn = '';
end
%this makes it PC and MAC compatible
if ~iscell(fn), fn = {fn}; end

for iFn = 1:length(fn)
    disp([9 'Importing from ', fn{iFn}, '... '])
    output(iFn) = importFromFile([pn, fn{iFn}]);
    disp([9 'Done'])
end


function output = importFromFile(fn)
%subfunction which does the actual job. For Axograph file header details
%see the Axograph User Manual
    fid = fopen(fn);

    output.filename = fn;
    OSType = fread(fid, 4, '*char','b')';
    if ~strcmp(OSType,'axgx')
        disp([9 'This type is not guaranteed to be supported... trying anyway!'])
    end
    fileFormat   = fread(fid, 1, 'int32','b')';
    nDatCol      = fread(fid, 1, 'int32','b')';

    for iYCol = 1:(nDatCol)
        nPoints  = fread(fid, 1, 'int32','b')';
        colType = fread(fid, 1, 'int32','b')';
        titlelen = fread(fid, 1, 'int32','b')';
        title = fread(fid, titlelen,'*char','b')';
        output.Channels(iYCol).title = title(2:2:length(title));
        
        switch colType
            case 4  % column type is short
                data(:, iYCol) = double(fread(fid,nPoints, 'int16','b'));
            case 5  % column type is long
                data(:, iYCol) = double(fread(fid,nPoints, 'int32','b'));
            case 6  % column type is float
                data(:, iYCol) = double(fread(fid,nPoints, 'float32','b'));
            case 7  % column type is double
                data(:, iYCol) = fread(fid,nPoints, 'double','b');
            case 9  % 'series' 
                 data0 = fread(fid,2, 'double','b');
                 data(:, iYCol) = (1:1:nPoints)*data0(2)+ data0(1);
            case 10 % 'scaled short'
                scale = fread(fid,1, 'double','b');
                offset = fread(fid,1, 'double','b');
                data0 = fread(fid,nPoints, 'int16','b');
                data(:, iYCol)  = double(data0)*scale + offset;
            otherwise
                disp(['Unknown column type:' num2str(output.Channels(iYCol).colType)  ' Cannot continue reading the file']);
                error('Exiting')
        end
        output.Channels(iYCol).data = data(:,iYCol);
    end

    fclose(fid);

