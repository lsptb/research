function MultiFunTA(Cfilenames,CfunHandle,nout,varargin)
% This function works on .axgd, .axgt and .mat files.
% Loads the requested files, performs the function(s) specified and returns
% 'Results', with the specified number of outputs per function to the
% workspace. The 'Results' can then be used as the input for this function,
% to perform further analyses or plotting (In such a case, Cfilenames =
% {'Results'}).
% Not familiar with function handles?  See below for a list of functions I
% use and example function handles to go with.
%
% INPUTS:
% Cfilenames - cell array of filenames, with extensions!
% CfunHandle - cell array of function handles
% nout       - number of output arguments for each function handle
%              (numerical array)
% varargin   - cell array containing any number of variables needed for the
%              specified functions. They must be the same length as the
%              file list, Cfilenames. i.e. They specify values or sets of
%              values for each file.
%
% OUTPUT:
% Results    - structure placed in the base workspace
%
%
%FUNCTION LIST AND EXAMPLE FUNCTION HANDLES (fh):
%
%O = PowerSpecTA(x,y,FreqRange,Bins,NormAbs,Notch)
%fh = @(x,y) PowerSpecTA(x,y,[10 80],8000,''Normalized'',[]);
%
%O = SpecGramTA(x,y,Smooth,Notch,Wind,WindOL,SmoothWindow)
%fh = @(x,y) SpecGramTA(x,y,''mtm'',[],8000,7800,2.5);
%
%O = XCorrTA(x,y1,y2,bpFiltParms,Notch,NormAbs)'
%fh = @(x,y1,y2) XCorrTA(x,y1,y2,[],[],''Absolute'');
%
%O = PhaseSyncTA(x,y1,y2,Wind,WindOL)
%fh = @(y1,y2) PhaseSyncTA(y1,y2,5000,500,490);
%
%O = IPSPsTA(x,yF,yIC,Size,Notch,method,FreqRange,Bins,SmoothWindow)
%fh = @(x,y1,y2) IPSPsTA(x,y1,y2,0.5,[],''pwelch'',[5 45],8000,[]);
%
%SpikeFieldCohTA(x,yF,yIC,bpFiltParms,Notch,varargin)
%fh = @(x,yF,yIC) SpikeFieldCohTA(x,yF,yIC,[10 80],''N'',section_start_sec,section_length_sec);
%
%XCorrGramTA(x,y1,y2,bpFiltParms,Notch,NormAbs)
%fh2 = @(x,y1,y2) XCorrGramTA(x,y1,y2,[10 80],[],''norm'');
%
%O = CharIhTA(x,y,offset_voltage,tonic_injected_current,sections_label,...
%sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
%fh = @(x,y,a,b,c,d,e,f,g) CharIhTA(x,y,a,b,c,d,e,f,g); (variables 'a'-'g'
%can be provided in a cell matrix and inputted as 'C{:}'.
%
%O = CharHyperpolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,...
%sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
%fh = @(x,y,a,b,c,d,e,f,g) CharHyperpolStepTA(x,y,a,b,c,d,e,f,g);
%(variables 'a'-'g' can be provided in a cell matrix and inputted as 'C{:}'.
%
%O = CharDepolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,...
%sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
%fh = @(x,yF,yIC,a,b,c,d) CharDepolTonicSpikesTA(x,yF,yIC,[],[],a,b,c,d)
%(variables 'a'-'d' can be provided in a cell matrix and inputted as 'C{:}'.
%
%O = CharDepolTonicSpikesTA(x,yF,yIC,bpFiltParms,Notch,offset_voltage,...
%tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec)
%fh = @(x,yF,yIC,a,b,c,d,e) CharDepolTonicSpikesTA(x,yF,yIC,[],[],a,b,c,d,e);
%(variables 'a'-'e' can be provided in a cell matrix and inputted as 'C{:}'.

%% get all data and input variables in place

% (for 'Ce', each cell is a different input for the function, with the
% length of each input equal to the number of files)


if ~iscell(CfunHandle)
    CfunHandle = {CfunHandle};
end

for kf = 1:length(CfunHandle)
    te = strsplit(func2str(CfunHandle{kf}),'(');
    te2 = cellfun(@(x) strsplit(x,')'),te,'Uni',0);
    funName{kf} = te2{2}{2};
    wf{kf} = which(funName{kf});
    if regexp(wf{kf},'built-in')
        error('You have entered a built-in function - these are inaccessible for use in this script. Sorry for any inconvenience caused.')
    end
end

%% EDITED BY JASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=Cfilenames;
[a,b,c]=fileparts(f{1});
if ~isempty(c)
  for k=1:length(f)
    [fp,Cfilenames{k},ext]=fileparts(f{k});
    Cfilenames_orig{k}=f{k};
  end
  try ext=strrep(ext,'.',''); end
else
  ext='var';
end
% e = cellfun(@(x) strsplit(x,'.'),Cfilenames,'Uni',0);
% if length(e{1})>1
%     ext = e{1}{end};
%     e2 = cellfun(@(x) x{1:end-1},e,'Uni',0);
%     e3 = cellfun(@(x) strjoin(x,''),e2,'Uni',0);
%     Cfilenames_orig = Cfilenames; Cfilenames = e3;
% else ext = 'var';
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ext
  
    case 'axgd'
        if isempty(Cfilenames_orig)
            O = importAxoX;  % Thank you Ed!!
        else O = importAxoX(Cfilenames_orig);
        end
        
        %add the extra variable list
        if nargin>3
            for kc = 2:2:length(varargin)
                %name = genvarname(['va' num2str(kc)]);
                name = genvarname([varargin{kc-1}]);
                for l = 1:length(O)
                    [O(l).Channels(end+1).title] = name;
                    if iscell(varargin{kc})
                        [O(l).Channels(end).data] = varargin{kc}{l};
                    else [O(l).Channels(end).data] = varargin{kc}(l,:);
                    end
                end
            end
        end
        
        for ko = 1:length(O)
            fn1 = strsplit(O(ko).filename,'\');
            fn2{ko} = strsplit(fn1{end},'.');
            for kc = 1:length(O(ko).Channels)
                terms = strrep(O(ko).Channels(kc).title,'\','');
                terms = strrep(terms,'(','');
                terms = strrep(terms,')','');
                varname{ko}{kc} = genvarname(terms);
                if ko==1, v{kc} = varname{ko}{kc}; end
                if strcmp(v{kc},varname{ko}{kc})
                    assignin('base','ko',ko);
                    assignin('base','temp',O(ko).Channels(kc).data);
                    evalin('base',[varname{ko}{kc} ' = temp;'])
                    for kf = 1:length(CfunHandle)
                        varval1{kf}{ko}{kc} = evalin('base','temp');
                    end
                else error('not all axograph files have the same name variables within. This will error when performing the requested functions.')
                end
            end
            clearvars -except Cfilenames CfunHandle cname O v varval1 fn2 nout ext wf funName
        end

    case 'mat'
        
        O = cellfun(@load,Cfilenames);
        
        %add the extra variable list
        if nargin>3
            for kc = 2:2:length(varargin)
                %name = genvarname(['va' num2str(kc)]);
                name = genvarname([varargin{kc-1}]);
                if iscell(varargin{kc})
                    [O(:).(name)] = deal(varargin{kc}{:});
                else[O(:).(name)] = deal(varargin{kc}(:));
                end
            end
        end
        
        for ko = 1:length(O)
            fn2{ko}{1} = genvarname(Cfilenames{ko});
            cname{ko} = fieldnames(O);
            for kc = 1:length(cname{ko})
                if ko==1, v{kc} = cname{ko}{kc}; end
                if strcmp(v{kc},cname{ko}{kc})
                    assignin('base','ko',ko);
                    assignin('base','temp',O(ko).(cname{ko}{kc}));
                    evalin('base',[cname{ko}{kc} ' = temp;'])
                    for kf = 1:length(CfunHandle)
                        varval1{kf}{ko}{kc} = evalin('base','temp');
                    end
                else error('not all files have the same name variables within. This will error when performing the requested functions.')
                end
            end
            clearvars -except Cfilenames CfunHandle cname O v varval1 fn2 nout ext wf funName
        end
        
    case 'axgt'
        
        O = cellfun(@importdata,Cfilenames_orig);
        for l = 1:length(O)
            rl = size(O(l).data,1); cl = ones(1,size(O(l).data,2));
            O(l).data = mat2cell(O(l).data,[rl],cl);
            for kc = 2:2:length(varargin)
                if iscell(varargin{kc})
                    O(l).data{end+1} = varargin{kc}{l};
                else O(l).data{end+1} = varargin{kc}(l,:);
                end
                %name = genvarname(['va' num2str(kc)]);
                name = genvarname([varargin{kc-1}]);
                O(l).textdata{end+1} = name;
            end
        end
        
        for ko = 1:length(O)
            fn2{ko}{1} = genvarname(Cfilenames{ko});
            for kc = 1:length(O(ko).textdata)
                cname{ko}{kc} = genvarname(regexprep(O(ko).textdata{kc},'[^\w'']',''));
            end
            for kc = 1:length(cname{ko})
                if ko==1, v{kc} = cname{ko}{kc}; end
                if strcmp(v{kc},cname{ko}{kc})
                    assignin('base','ko',ko);
                    assignin('base','temp',O(ko).data{kc});
                    evalin('base',[cname{ko}{kc} ' = temp;'])
                    for kf = 1:length(CfunHandle)
                        varval1{kf}{ko}{kc} = evalin('base','temp');
                    end
                else error('not all files have the same name variables within. This will error when performing the requested functions.')
                end
            end
            clearvars -except Cfilenames cname CfunHandle O v varval1 fn2 nout ext wf funName
        end
        
    case 'var'
        
        name1 = evalin('base',['fieldnames(' Cfilenames{1} ');']);
        name2 = evalin('base',['fieldnames(' Cfilenames{1} '.' name1{1} ');']);
        for k = 1:length(name2)
            O(k) = evalin('base',[Cfilenames{1} '.' name1{1} '.' name2{k} ';']);
        end
        Cfilenames = name2;
        for ko = 1:length(O)
            fn2{ko}{1} = genvarname(Cfilenames{ko});
            cname{ko} = fieldnames(O{1});
            for kc = 1:length(cname{ko})
                if ko==1, v{kc} = cname{ko}{kc}; end
                if strcmp(v{kc},cname{ko}{kc})
                    assignin('base','ko',ko);
                    assignin('base','temp',O{ko}.(cname{ko}{kc}));
                    evalin('base',[cname{ko}{kc} ' = temp;'])
                    for kf = 1:length(CfunHandle)
                        varval1{kf}{ko}{kc} = evalin('base','temp');
                    end
                else error('not all files have the same name variables within. This will error when performing the requested functions.')
                end
            end
            clearvars -except Cfilenames cname CfunHandle O v varval1 fn2 nout ext wf funName
        end
        
end

%% request input labels

inpts = cellfun(@(x) num2cell(zeros(1,nargin(x))),CfunHandle,'Uni',0);
for k = 1:length(inpts)
    [inpts{k}{:}] = deal(''); inpts{k} = inpts{k}(:);
    hf(k) = figure('Position',[100 100 500 600]);
%     if strcmp('var',ext)
%         t1 = strsplit(wf{k},'(');
%         t2 = strsplit(t1{end},')');
%         t3 = strcat(t2(1),{'.m'}); wf{k} = t3{1};
%     end
    FC = strsplit(fileread(wf{k}),')'); FC(2:end) = []; FC{1} = strcat(FC{1},')');
    vterms = cellfun(@(x) strsplit(x,'_'),v,'Uni',0);
    for iv = 1:length(vterms), v{iv} = strjoin(vterms{iv},''); end
    annotation(hf(k),'textbox',[0.06 0.6 0.9 0.3],...
        'String',strcat([FC ' ' v]),'FitBoxToText','off','LineStyle','none','FontSize',8);
    uitable('Parent',hf(k),'ColumnEditable',true,'ColumnFormat',{'char'},... 
        'ColumnName',{[func2str(CfunHandle{k}) ' Inputs:']},'Data',inpts{k},...
        'Position',[20 20 450 300],'CellEditCallback',@uicall);
    uiwait(hf(k))
end

%% perform analysis

for kf = 1:length(CfunHandle)
    disp(['Performing Function ' num2str(kf) '..'])
    hwb = waitbar(0,'');
    for ko = 1:length(O)
        tic
        %try
        RCell.(funName{kf}).(genvarname(fn2{ko}{1})) = num2cell(zeros(1,nout(kf)));
            [RCell.(funName{kf}).(genvarname(fn2{ko}{1})){:}] = feval(CfunHandle{kf},varval2{kf}{ko}{:});
            for ka = 1:length(RCell.(funName{kf}).(genvarname(fn2{ko}{1})))
                name = genvarname(['output_' num2str(ka)]);
                Results.(funName{kf}).(genvarname(fn2{ko}{1})).(name) = RCell.(funName{kf}).(genvarname(fn2{ko}{1})){ka};
            end
            waitbar(ko/length(O),hwb,['Completed ''' Cfilenames{ko} ''''])
        %catch %#ok<CTCH>
        %    disp(['Failed analysis of ''' Cfilenames{ko} ''''])
        %end
        toc
    end
    close(hwb)
end
if sum(nout(kf)>0)
    assignin('base','Results',Results);
    explorestructTA(Results)
end
disp('Analysis complete')
evalin('base',['clear ','ko'])
evalin('base',['clear ','temp'])
if exist('cname','var')
    for koo = 1:length(cname)
        for kcc = 1:length(cname{koo})
            evalin('base',['clear ',cname{koo}{kcc}])
        end
    end
end

%% nested functions
    function uicall(~,e)
        inpts{k}{e.Indices(1),e.Indices(2)} = e.NewData;
        li = strcmp(e.NewData,v);
        if sum(li)<1
            for ko2 = 1:length(O)
                val = str2double(e.NewData);
                if isempty(val)
                    valstr = e.NewData;
                    varval2{k}{ko2}{e.Indices(1),e.Indices(2)} = valstr;
                else varval2{k}{ko2}{e.Indices(1),e.Indices(2)} = val;
                end
            end
        else
            for ko2 = 1:length(O)
                if strcmp('var',ext)
                    varval2{k}{ko2}{e.Indices(1),e.Indices(2)} = varval1{k}{ko2}{li}{1};
                else varval2{k}{ko2}{e.Indices(1),e.Indices(2)} = varval1{k}{ko2}{li};
                end
            end
        end
    end


end

