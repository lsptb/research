codepath='/project/crc-nak/sherfey/code/research/external/tallie';
addpath(genpath(codepath));

%% load experimental data (use Tallie's code)
% ...

% get parameters for cell characterization
parms = getCharParms(parmfile,'exp');

% characterize cells
funnames = {'CharHyperpolStepTA','CharDepolStepTA','CharDepolTonicSpikesTA'};
res = {};
for f=1:length(funnames)
  args = struct2cell(parms.(funames{f}));
  res{f} = feval(funnames{f},args{:});
end

% construct distributions
% ...


%% load simulated data and model specification

% this function should be called by a higher level function that takes a
% list of files from a simulation batch and processes them all, saving
% post-processed results into a batch analysis directory.

datafile = '/project/crc-nak/sherfey/code/tests/cellchar/20140215-100816/E-multiplicity/data/E-multiplicity_1_sim_data.mat';
load(datafile,'sim_data','spec','parms');
cellname = 'E'; thiscell = 1; datavar = [cellname '_V'];

cellind = find(strcmp(cellname,{spec.entities.label}));
dataind = find(strcmp(datavar,{sim_data(cellind).sensor_info.label}));

t = sim_data(cellind).epochs.time;
intra = sim_data(cellind).epochs.data(dataind,:,thiscell);
field = mean(sim_data(cellind).epochs.data(dataind,:,:),3);

% get parameters for cell characterization
aparms = getCharParms(spec,'sim','cellind',cellind,'mechanism','iStepProtocol');

% characterize cells
funnames = {'CharHyperpolStepTA','CharDepolStepTA','CharDepolTonicSpikesTA'};
dataargs = {{t,intra}           ,{t,intra}        ,{t,field,intra}};
res = {};
for f=1:length(funnames)
  args = struct2cell(aparms.(funnames{f}));
  args = {dataargs{f}{:},args{:}};
  res{f} = feval(funnames{f},args{:});
end

% ...


