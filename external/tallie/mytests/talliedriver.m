codepath='/project/crc-nak/sherfey/code/research/external/tallie';
datapath='/project/crc-nak/sherfey/projects/rhythms/rat/tallietoolbox';
filefile=fullfile(datapath,'/Files & Function Handle.mat');
parmfile=fullfile(datapath,'egparms.xlsx');

addpath(genpath(codepath));
cd(datapath);

parms = GetParmsTA(parmfile);
%{
    What field do you want to select by? 'Project'
    What value should this field have? 'AC Characterization'
    How many fields do you want to collate into easier-to-access matrices? 4
    Fieldname 1: 'tonic_injected_current'
    'num' or 'str': 'num'
    Fieldname 2: 'sections_label'
    'num' or 'str': 'num'
    Fieldname 3: 'sections_start_sec'
    'num' or 'str': 'num'
    Fieldname 4: 'sections_length_sec'
    'num' or 'str': 'num'
    'The requested variables are held in the workspace, both as'
    'individual variables and as a cell array of variables, 'Vars''
    'and their labels, 'Vars_labels'. These are concatenated into 'VarsC'.'
%}

fileinfo=load(filefile);
files=var.FilesC; for k=1:length(files),files{k}=fullfile(datapath,files{k}); end
fh=var.fh;

MultiFunTA(files',fh,1,VarsC{:}) % (Cfilenames,CfunHandle,nout,varargin)
  % note: had to modify the filename parsing b/c of an error (see lines 88-109)
  
%{
Tallie: then a box comes up listing the variable available for the function.  
You want 'TimeS' for 'x' and 'AnalogInput0V' for 'yIC', and then the last 
three listed variables for the three remaining spaces in the table.

fn: @(x,yIC,tc,sl,sss,sls)PowerSpecWobbleTA(x,yIC,[2,45],[],3.5,tc,sl,sss,sls)
parms to enter into uitable:
  1 = TimeS
  2 = AnalogInput0V
  3 = tonicinjectedcurrent
  4 = sectionslabel
  5 = sectionsstartsec
  6 = sectionslengthsec
note: then close the window to continue

After running my function you'll get my adapted vertsion of explorstruct open up.  
You can export individual variables, or alteratively, all of a chosen variable 
for every file by filling in your function handle name in the box that comes up.

All results are saved in a single structure 'Results':

Results.PowerSpecWobbleTA.
    ci56121112010: [1x1 struct]
    ci57121112023: [1x1 struct]
    ci77121210005: [1x1 struct]
    ci81121211022: [1x1 struct]
Results.PowerSpecWobbleTA.ci56121112010.output_1.
                          raw: [1908763x1 double]
       tonic_injected_current: NaN
                    FreqRange: [2 45]
               sections_label: [0 100 200 300 400 500 600 700 800 900 100]
           sections_start_sec: [6.0000e-04 25 45 65 85 125 163 200 240 285 322]
          sections_length_sec: [15 15 15 15 15 15 17 20 20 15 18]
                data_sections: [1x1 struct]
                    AreaPower: [1x1 struct]
                    PeakPower: [1x1 struct]
                      OscFreq: [1x1 struct]
                          Pxx: [1x1 struct]
                            f: [1x1 struct]
    total_length_sections_sec: [15 15 15 15 15 14.9998 17 9.5498 2.7122]
                     Pxx_mean: [1x1 struct]
                       f_mean: [1x1 struct]
               AreaPower_mean: [1x1 struct]
               PeakPower_mean: [1x1 struct]
                 OscFreq_mean: [1x1 struct]
                         data: [1x1 struct]
                            t: [1x1 struct]
                           Fs: 5000
                 Pxx_HzPerBin: 0.2500
                LessThan2Secs: 'True'
%}


%{
  
NEXT:
  1. test plot funcs (MultiFigTA.m)
  2. test other analysis funcs (eg, PowerSpecTA.m, XCorrGramTA.m)
      - see /project/crc-nak/sherfey/code/research/external/tallie/funTA
      - also see Tallie comments under specific analyses below.
  3. understand how parameter setting/variation works
  4. figure out how to plug in sim data or extract analysis/plot guts for own funcs
  5. incorporate into cellmodeler & netmodeler
  6. then incorporate simstudy controls and apply get_search_space to form sets
  of simulations as well as sets of analyses (w/ cluster handling + qmatjobs 
  scripts copied from mmil, and structured prefixes/rootoutdirs like Shane)
  
%}

%% 1. test plot funcs
% MultiFigTA(FuncHandle,Vars)
  % FuncHandle = @plot
  % VarsStruct.f1...fn = {[input1] [input2]}
  % VarsCell{{1x2}...{1x2}} = {input1} {input2}
  % VarsDouble[] = columns to subplot individually

% ex) cell array of data
figure
h = @plot;
s = Results.PowerSpecWobbleTA.ci56121112010.output_1;
MultiFigTA(h,{s.raw, s.data.sect_100, s.data.sect_200});
% had to comment out 'fm'. is this an external function i don't have???

% ex) structure with data fields
figure
fields = {'sect_100','sect_200','sect_300'};
s2 = rmfield(s.data,setdiff(fieldnames(s.data),fields));
MultiFigTA(h,s2);

% ex) data matrix
figure
M = [s2.sect_100 s2.sect_200 s2.sect_300];
MultiFigTA(h,M);


%% 2. test other analysis funcs (eg, PowerSpecTA.m, XCorrGramTA.m)
% - see /project/crc-nak/sherfey/code/research/external/tallie/funTA
%{
Comments on analyses of interest --
Tallie: I use power spectra and cross correlation functions with working with LFPs. 
For intracels, I use all the cell charaterization protocols at the minute 
as I'm trying to decipher a meaningful selection of parameters I might use 
to automate the clustering cell types. I use the power spectra wobble 
function for looking at STOs. I will use IPSP and EPSP programmes (when I 
get to it), along with cross correlations and phase/sync functions to look 
at the relationship of the cell inputs to the field.  Finally, I am also 
using CSD analysis with multielectrode data using CSDplotter (downloadable).
...I didn't mention the spike/field coherence stuff, but I'll use that too 
as an assessment of spike-time dependence on drive to the cell.

LFP: PowerSpecTA, XCorrTA, XCorrGramTA
Vm: CharHyperpolStepTA, CharDepolTonicSpikesTA, CharDepolTonicSpikesTA
STOs: PowerSpecWobbleTA
Synapse/Field relationship: PhaseSyncTA, IPSPsTA, XCorrTA
Spike/Field relationship: SpikeFieldCohTA
CSD: CSDplotter

%}

% LFP:

%O = PowerSpecTA(x,y,FreqRange,Bins,NormAbs,Notch)
tic; figure; fh = @(x,y) PowerSpecTA(x,y,[10 80],8000,'Normalized',[]);
MultiFunTA(files',fh,1,VarsC{:}); toc
  o=Results.PowerSpecTA.ci56121112010.output_1;
  MultiFigTA(@plot,);
  figure; plot(o.f,log10(o.Pxx)); xlim([0 100]);

%O = SpecGramTA(x,y,Smooth,Notch,Wind,WindOL,SmoothWindow)
tic; figure; fh = @(x,y) SpecGramTA(x,y,'mtm',[],8000,7800,2.5);
MultiFunTA(files',fh,1,VarsC{:}); toc
  o=Results.SpecGramTA.ci56121112010.output_1; specgram=Results;
  y=o.yo; % freq x time
  yn=(y-repmat(mean(y,2),[1 size(y,2)]))./repmat(std(y,[],2),[1 size(y,2)]);
  figure
  subplot(3,1,1); imagesc(o.to,o.fo,yn); axis xy; title('spectrogram (z-score)');
  colormap(flipud(colormap('gray'))); colorbar;
  subplot(3,1,2); imagesc(o.to,o.fo,yn); axis xy; ylim([0 100]); colorbar;   
  subplot(3,1,3); imagesc(o.to,o.fo,yn); axis xy; colorbar;
  ylim([0 100]); xlim([.8 .9]*o.to(end));
  xlabel('time (s)'); ylabel('freq (Hz)'); title('spectrogram (z-score)');
  
% Vm:

%O = CharHyperpolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,...
%sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
figure; tic; fh = @(x,y,a,b,c,d,e,f,g) CharHyperpolStepTA(x,y,a,b,c,d,e,f,g);
MultiFunTA(files',fh,1,VarsC{:}); toc
  %   TimeS
  %   AnalogInput0V
  %   NaN [Offset]
  %   tonicinjectedcurrent
  %   sectionslabel
  %   sectionsstartsec
  %   sectionslengthsec
  %   NaN [baseline_start_sec]
  %   NaN [baseline_length_sec]

%O = CharDepolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,...
%sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
figure; tic; fh = @(x,yF,yIC,a,b,c,d) CharDepolTonicSpikesTA(x,yF,yIC,[],[],a,b,c,d)
MultiFunTA(files',fh,1,VarsC{:}); toc

%O = CharDepolTonicSpikesTA(x,yF,yIC,bpFiltParms,Notch,offset_voltage,...
%tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec)
figure; tic; fh = @(x,yF,yIC,a,b,c,d,e) CharDepolTonicSpikesTA(x,yF,yIC,[],[],a,b,c,d,e);
MultiFunTA(files',fh,1,VarsC{:}); toc
  
% Synapse/Field relationship:
  
%O = PhaseSyncTA(x,y1,y2,Wind,WindOL)
figure; tic; fh = @(y1,y2) PhaseSyncTA(y1,y2,5000,500,490);
MultiFunTA(files',fh,1,VarsC{:}); toc

%O = IPSPsTA(x,yF,yIC,Size,Notch,method,FreqRange,Bins,SmoothWindow)
figure; tic; fh = @(x,y1,y2) IPSPsTA(x,y1,y2,0.5,[],'pwelch',[5 45],8000,[]);
MultiFunTA(files',fh,1,VarsC{:}); toc

% Spike/Field relationship:

%SpikeFieldCohTA(x,yF,yIC,bpFiltParms,Notch,varargin)
figure; tic; fh = @(x,yF,yIC) SpikeFieldCohTA(x,yF,yIC,[10 80],'N',section_start_sec,section_length_sec);
MultiFunTA(files',fh,1,VarsC{:}); toc

% STOs: PowerSpecWobbleTA

% Ih?:
%O = CharIhTA(x,y,offset_voltage,tonic_injected_current,sections_label,...
%sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
figure; tic; fh = @(x,y,a,b,c,d,e,f,g) CharIhTA(x,y,a,b,c,d,e,f,g); 
MultiFunTA(files',fh,1,VarsC{:}); toc
% (variables 'a'-'g' can be provided in a cell matrix and inputted as 'C{:}'.

% FAILED:
% %XCorrGramTA(x,y1,y2,bpFiltParms,Notch,NormAbs)
% figure; tic; fh = @(x,y1,y2) XCorrGramTA(x,y1,y2,[10 80],[],'norm');
% MultiFunTA(files',fh,0,VarsC{:}); toc
% 
% %O = XCorrTA(x,y1,y2,bpFiltParms,Notch,NormAbs)'
% figure; tic; fh = @(x,y1,y2) XCorrTA(x,y1,y2,[],[],'Absolute'); 
% MultiFunTA(files',fh,1,VarsC{:}); toc

% all functions togetherfh={};
fh{end+1} = @(x,y) PowerSpecTA(x,y,[10 80],8000,'Normalized',[]);
fh{end+1} = @(x,y) SpecGramTA(x,y,'mtm',[],8000,7800,2.5);
% fh{end+1} = @(x,y1,y2) XCorrTA(x,y1,y2,[],[],'Absolute'); 
fh{end+1} = @(y1,y2) PhaseSyncTA(y1,y2,5000,500,490);
fh{end+1} = @(x,y1,y2) IPSPsTA(x,y1,y2,0.5,[],'pwelch',[5 45],8000,[]);
fh{end+1} = @(x,yF,yIC) SpikeFieldCohTA(x,yF,yIC,[10 80],'N',section_start_sec,section_length_sec);
fh{end+1} = @(x,y,a,b,c,d,e,f,g) CharIhTA(x,y,a,b,c,d,e,f,g); 
fh{end+1} = @(x,y,a,b,c,d,e,f,g) CharHyperpolStepTA(x,y,a,b,c,d,e,f,g);
fh{end+1} = @(x,yF,yIC,a,b,c,d) CharDepolTonicSpikesTA(x,yF,yIC,[],[],a,b,c,d)
fh{end+1} = @(x,yF,yIC,a,b,c,d,e) CharDepolTonicSpikesTA(x,yF,yIC,[],[],a,b,c,d,e);
MultiFunTA(files',fh,ones(1,length(fh)),VarsC{:}); toc

%% figure out how to plug in sim data or extract analysis/plot guts for own funcs
O = importAxoX(files);

% manual visualization
nd=length(O);
nc=cellfun(@length,{O.Channels});
nt=arrayfun(@(x)length(x.Channels(1).data),O);
nsamp=max(nt);
data=zeros(sum(nc),nsamp);
cnt=0;
for i=1:nd
  for j=1:nc(i)
    cnt=cnt+1;
    dat=O(i).Channels(j).data;
    data(cnt,1:length(dat))=dat;
  end
end
data=ts_matrix2data(data,'sfreq',1/(data(1,2)-data(1,1)));
visualizer(data); % check chan2 @ 330-340 (clear ~8-10Hz STO; img view)

% manual spectral analysis
fh = @(x,y) PowerSpecTA(x,y,[10 80],8000,'Normalized',[]);
clear res
figure; 
for i=1:nd
  t=O(i).Channels(1).data; 
  V=O(i).Channels(2).data; 
  I=O(i).Channels(3).data; 
  res(i) = feval(fh,t,V);
  subplot(3,nd,i); plot(res(i).f,log10(res(i).Pxx)); xlim([0 100]);
  subplot(3,nd,i+nd); plot(t,V);
  subplot(3,nd,i+2*nd); plot(t,I);
end

% given simdata:
pop=1; var=1;
t=simdata(pop).epochs.time;
V=simdata(pop).epochs.data(var,:);
res = feval(fh,t,V);
  
%% test other data sets (other formats; eg, .axgt, .mat)
file = {'/project/crc-nak/sherfey/projects/rhythms/rat/predelta 020.txt.axgt'};
% note: this is data I collected during my first visit to Newcastle
fh = @(x,y) PowerSpecTA(x,y,[10 80],8000,'Normalized',[]);

% MultiFunTA load and analysis
MultiFunTA(file,fh,1); % (Cfilenames,CfunHandle,nout,varargin)
o=Results.PowerSpecTA.predelta0200x2Etxt.output_1;
figure; plot(o.f,log10(o.Pxx)); xlim([0 100]);  

% manual load and analysis
axgt = cellfun(@importdata,file);
data = axgt.data';
data=ts_matrix2data(data,'time',data(1,:));
visualizer(data);

t=axgt.data(:,1);
V=axgt.data(:,2);
I=axgt.data(:,3);
fh = @(x,y) PowerSpecTA(x,y,[10 80],8000,'Normalized',[]);
res = feval(fh,t,V);
figure
subplot(3,1,1); plot(res.f,log10(res.Pxx)); xlim([0 100]);
subplot(3,1,2); plot(t,V);
subplot(3,1,3); plot(t,I);

% 3 approaches to plugging simdata into MultiFunTA():
% - modify MultiFunTA to load simdata and organize into proper structure
% - save simdata in Tallie format for use w/ MultiFunTA "as is"
% - skip MultiFunTA and pass simdata (t,x) to the desired fh directly

%%

codepath='/project/crc-nak/sherfey/code/research/external/tallie';
datapath='/project/crc-nak/sherfey/projects/rhythms/rat/cell-characterization/hyper and depol steps/hyper/';
parmfile='/project/crc-nak/sherfey/projects/rhythms/rat/cell-characterization/Parms hyperpol steps blockers temp.xlsx';
addpath(genpath(codepath));
cd(datapath);

parms = GetParmsTA(parmfile);
%{
What field do you want to select by? 'trace_description'
What value should this field have? 'hyperpol steps'
How many fields do you want to collate into easier-to-access matrices? 7
Fieldname 1: 'offset_potential'
'num' or 'str': 'num'
Fieldname 2: 'tonic_injected_current'
'num' or 'str': 'num'
Fieldname 3: 'sections_label'
'num' or 'str': 'num'
Fieldname 4: 'sections_start_sec'
'num' or 'str': 'num'
Fieldname 5: 'sections_length_sec'
'num' or 'str': 'num'
Fieldname 6: 'baseline_start_sec'
'num' or 'str': 'num'
Fieldname 7: 'baseline_length_sec'
'num' or 'str': 'num'
    'Done.'
%}
%fh = CharHyperpolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec);

fh=@(x,y,a,b,c,d,e,f,g)CharHyperpolStepTA(x,y,a,b,c,d,e,f,g);
f=dir; 
f={f(~[f.isdir]).name};
f=f(~cellfun(@isempty,regexp(f,'ci.*.axgt')));

sel=1:3;
fs=f(sel);
vs=VarsC; vs(2:2:end)=cellfun(@(x)x(sel),vs(2:2:end),'uni',0);
MultiFunTA(fs,fh,1,vs{:}) % (Cfilenames,CfunHandle,nout,varargin)

% output_1.step_sections (@plot)
% export multiples: Results.CharHyperpolStepTA
figure; MultiFigTA(@plot,V);
% output_1.Ih_Peak_mV (@bar)
% export multiples: Results.CharHyperpolStepTA
figure; MultiFigTA(@bar,V);


