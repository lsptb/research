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


%% Specific analyses
%{
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
%}




