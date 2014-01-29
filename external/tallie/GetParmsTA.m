function Parms = GetParmsTA(xlsFile)
% Reads the xls spreadsheet of parameter values (the first column being the
% headers and columns 8 onwards being the parameter values for each
% reading) and organizes that data into a structure.  This structure is
% should be all of the metadata necessary for any analyses on a particular
% set of traces. (Columns 2-7 are not included as they are example
% parameter values for some standard readings.)
%
% It then asks you which files you want to work with (selected out by a
% specified field contents), and whether you want certain variables pulling
% out for easier access. If you do, they are created in the base workspace.
%
% INPUTS:
% xlsFile   - 'Parameters.xlsx'
%

Parms = ParamXls1TA(xlsFile);
fieldnames(Parms)

fieldname = input('What field do you want to select by? ');
val = input('What value should this field have? ');
Parms = SelectParms2TA(Parms,fieldname,val); assignin('base','Parms',Parms)

kn = input('How many fields do you want to collate into easier-to-access matrices? ');
assignin('base','VarsCell',{})
Vars_labels = {};
for k = 1:kn
    fieldname = input(['Fieldname ' num2str(k) ': ']);
    assignin('base','fieldname',fieldname)
    NumStr = input('''num'' or ''str'': ');
    assignin('base','NumStr',NumStr)
    Parms = evalin('base','Parms');
    fnv = CollateParms3TA(Parms,fieldname,NumStr)';
    assignin('base',fieldname,fnv)
    Vars = evalin('base',['cat(1,VarsCell,{' fieldname '});']);
    Vars_labels = cat(1,Vars_labels,{fieldname});
    assignin('base','VarsCell',Vars);
end
VC = {};
for k = 1:length(Vars_labels)
    VC{end+1} = cat(2,Vars_labels(k),Vars(k));
end
VarsC = {};
for k = 1:length(VC)
    VarsC = cat(2,VarsC,VC{k});
end
assignin('base','Vars',Vars)
assignin('base','Vars_labels',Vars_labels)
assignin('base','VarsC',VarsC)
evalin('base','clear fieldname NumStr VarsCell')
disp({'Done.';...
    'The requested variables are held in the workspace, both as';...
    'individual variables and as a cell array of variables, ''Vars''';...
    'and their labels, ''Vars_labels''. These are concatenated into ''VarsC''.';})

end

