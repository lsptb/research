function Parms = ParamXls1TA(xlsFile)
% Reads the xls spreadsheet of parameter values (the first column being the
% headers and columns 8 onwards being the parameter values for each
% reading) and organizes that data into a structure.  This structure is
% should be all of the metadata necessary for any analyses on a particular
% set of traces. (Column 2 is not included as it is a dummy example.)
%
% INPUTS:
% xlsFile   - 'Parameters.xlsx'
%

% read Parameters.xlsx and form a structure
[~,~,r] = xlsread(xlsFile);
Parms = cell2struct(r(:,8:end),r(:,1),1);

% field names for any data with a comma delimiter are split into arrays
sf = fieldnames(Parms);
for i = 1:length(Parms)
    for k = 1:length(sf)
        if ~isnumeric(Parms(i).(sf{k}))
            temp = Parms(i).(sf{k});
            ch = strsplit(temp,',');
            if length(ch)>1
                Parms(i).(sf{k}) = ch;
            end
        end
        clear temp ch
    end
end

end

