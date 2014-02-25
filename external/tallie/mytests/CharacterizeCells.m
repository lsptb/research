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

% construct feature distributions
% ...

% calculate mean and standard deviations for z-scoring simulated responses
% ...


%% load simulated data and model specification

% this function should be called by a higher level function that takes a
% list of files from a simulation batch and processes them all, saving
% post-processed results into a batch analysis directory.

%datafile = '/project/crc-nak/sherfey/projects/rhythms/models/tests/20140223-163742/E-multiplicity/data/E-multiplicity_1_sim_data.mat';
%datafile = '/project/crc-nak/sherfey/projects/rhythms/models/tests/20140223-163742/E-multiplicity/data/E-multiplicity_2_sim_data.mat';
datafile = '/project/crc-nak/sherfey/projects/rhythms/models/tests/20140223-165742/E2-multiplicity_E-multiplicity/data/E2-multiplicity2_E-multiplicity2_sim_data.mat';

load(datafile,'sim_data','spec','parms');
var='V'; plot_flag=0;

tic
res = {};
cellnames={spec.entities.label};
for thispop = 1:length(cellnames)

  cellname = cellnames{thispop}; %'E'; 
  datavar = [cellname '_' var]; % 'E_V';

  cellind = find(strcmp(cellname,{spec.entities.label}));
  dataind = find(strcmp(datavar,{sim_data(cellind).sensor_info.label}));

  for thiscell = 1:spec.entities(cellind).multiplicity

    t = sim_data(cellind).epochs.time;
    intra = double(sim_data(cellind).epochs.data(dataind,:,thiscell)/1000);
    field = double(mean(sim_data(cellind).epochs.data(dataind,:,:),3)/1000);
    sz=200/1000; % artifact size, [V]
    prestep_time_offset = -.01; % shift sections and artifacts so onset artifact can be detected by tallie's cell characterization scripts

    % get parameters for cell characterization
    aparms = getCharParms(spec,'sim','cellind',cellind,'mechanism','iStepProtocol','prestep_time_offset',prestep_time_offset);
    funnames = {'CharHyperpolStepTA','CharDepolStepTA','CharDepolTonicSpikesTA'};

    % add artifact to hyperpolarization and depolarization steps
    for f=1:2
      % make hyperpol artifacts up=>down; depol down=>up
      p = spec.entities(cellind).mechs(1).params;
      tstart=aparms.(funnames{f}).sections_start_sec - prestep_time_offset;
      tinf=tstart+((0:p.nsteps*p.nsections-1)*p.isi)/1000;
      tind=round(tinf/(t(2)-t(1)));
      intra(tind)=sz*((-1)^(f-1)); % onset artifact
      tinf=tinf+p.steptime/1000;
      tind=round(tinf/(t(2)-t(1)));
      intra(tind)=-sz*((-1)^(f-1)); % offset artifact
      %figure('position',[25 420 1550 250]); plot(t,intra);
    end

    % characterize cells
    dataargs = {{t,intra}           ,{t,intra}        ,{t,field,intra}};
    for f=1:length(funnames)
      args = struct2cell(aparms.(funnames{f}));
      args = {dataargs{f}{:},args{:}};
      res{thiscell,f,thispop} = feval(funnames{f},args{:},plot_flag);
    end % end loop over cell characterization functions

  end % end loop over cells
  
end % end loop over populations
toc

%% collect response features
hyper=1; depol=2; tonic=3; thispop=1; thiscell=1;

p = spec.entities(1).mechs(1).params;
Fs = sim_data(thispop).sfreq;
tpulse = 0:1/Fs:(p.steptime/1000); %tstep = (0:size(hypdat,1)-1)/Fs;

% construct drive vector (current input per step)
amp = p.stepsize;
nsect = p.nsections;
nstep = p.nsteps;
%inputs = zeros(nstep,2*nsect); % step x sect(hyper,depol)
%for k=1:nsect, inputs(:,k) = (k*amp)*ones(1,nstep); end
%for k=(nsect+1):(2*nsect), inputs(:,k) = -(k*amp)*ones(1,nstep); end

o1 = res{thiscell,hyper,thispop};
o2 = res{thiscell,depol,thispop};
o3 = res{thiscell,tonic,thispop};

% passive membrane properties (Rin, gleak, Vrest, Cm)
if ~isnumeric(o1.Resistance_Mohms), Rin1=[]; else Rin1=o1.Resistance_Mohms; end
if ~isnumeric(o2.Resistance_Mohms), Rin2=[]; else Rin2=o2.Resistance_Mohms; end
Vrest1 = o1.Baseline_mV;
Vrest2 = o2.Baseline_mV;
% gL=1/Rin? Cm?

% calculate I/V and f/I curves from hyperpol and depol steps
% I/V
hypdat = o1.step_sections; % time x sect, mean step per section
depdat = o2.y_NoSpike_sect; % time x step x sect, step per section
sel = round([.2 .8]*length(tpulse));
sel = sel(1):sel(2);
Vlo = mean(hypdat(sel,:),1); Vlo(isnan(Vlo))=[];
Ilo = -amp*(1:length(Vlo));
Vhi = mean(depdat(sel,:),1); Vhi(isnan(Vhi))=[];
Ihi = amp*(1:length(Vhi));
V = [Vlo Vhi]; IVv=V;
I = [Ilo Ihi]; IVi=I;
P = polyfit(I,V,1);
IVslope = P(1);
IVinter = P(2);
figure; plot(I,V,'b*--','markersize',10); xlabel('current'); ylabel('voltage');

% f/I
rates = 1./cellfun(@median,o2.Spikes_ISI_median);
Irate = (1:nsect)*amp;
rateIC = nan(size(rates)); % instantaneous rate at second spike
rateSS = nan(size(rates)); % steady state firing rate
for i=1:nsect
  if ~isnan(rates(i))
    spks = o2.Spikes_InstFreq{i}; % o2.Spikes_InstFreq{sect}{step}
    rateIC(i) = median(cellfun(@(x)x(1),spks));
    rateSS(i) = median(cellfun(@(x)x(end-1),spks));
  end
end
P = polyfit(Irate,rates,1);
FIslope = P(1);
FIinter = P(2);

% h-channel kinetic variability
% ...

% ISI histogram / spike accommodation (minimize MSE)
% ...

% collapse response properties into scalar features
% ...

% define response feature vector
features = [Rin1 Rin2 Vrest1 Vrest2 IVslope IVinter FIslope FIinter];
feature_labels = {'Rin1','Rin2','Vrest1','Vrest2','IVslope','IVinter','FIslope','FIinter'};
feature_weights = ones(size(features));

% load or calculate experimental mean and standard deviations per feature
% note: this will probably be accomplished at the beginning of this script.
% for now set to 0s and 1s.
expt_mu = zeros(size(features));
expt_sd = ones(size(features));

% calculate z-scores on simulated feature vectors
zfeatures = (features - expt_mu) ./ expt_sd;

% calculate mean square feature z-scores (MSFZ)
MSFZ = mean(zfeatures.^2);
MSWFZ = mean((zfeatures.*feature_weights).^2); % weighted version

% save response features to batch analysis directory
% save(featurefile,'features','feature_labels','expt_mu','expt_sd','zfeatures','MSFZ','feature_weights','MSWFZ','IVi','IVv');%,'I','V'
% ...

% This function's caller should then load all response features
% across all simulations, construct the joint parameter space and calc
% fitness functions over it. lastly it should analyze the fitness surface
% and do model selection...  
%   e.g., 1) fitness = MSFZ - f(#params)
%         2) select models with fitness in lower 5%-ile
%         3) inspect mechanisms and parameter manifolds of selected models

%%
%{
cell model constraints:
- passive membrane properties (Rin, gleak, Vrest, Cm)
- calculate I/V and f/I curves from hyperpol and depol steps
- h-channel kinetic variability
- ISI histogram / spike accommodation (minimize MSE)

Cells:
passive membrane parameters:
Rin = CharHyperpolStepTA().Resistance_Mohms
Rin = CharDepolStepTA().Resistance_Mohms
gleak = 1/Rin
Vrest =CharHyperpolStepTA().Baseline_mV
Vrest = CharDepolStepTA().Baseline_mV
Cm?

I/V curve:
foreach section (I): V(I) = avg(V| end interval)
     CharHyperpolStepTA().step_sections
     CharDepolStepTA().step_sections
construct complete I/V curve by combining hyperpol and depol curves.

f/I curve:
CharDepolStepTA().Spikes_ISIs

h-channel:
CharHyperpolStepTA().IsThereIhYNeach -- proportion of cells with ih
CharHyperpolStepTA().output.IhDecayToThird_ms
CharHyperpolStepTA().output.Ih_mV
use their distributions to set population heterogeneity of gh and gtau

ISI distribution + spike accommodation (aka: frequency adaptation):
CharDepolStepTA().Spikes_ISIs -- frequency adaptation
CharDepolTonicSpikesTA().ISIs -- distribution

Spikes
shape -- CharDepolTonicSpikesTA().Spike_ave
SpikeAHP_HalfWidth_ms - AHP width at half-height (of AHP)
Spike_Amplitude - look at amplitude distribution
Spike_Width
rate -- CharDepolTonicSpikesTA().SpikeRate_persec

STOs
PowerSpecWobbleTA()

Networks:
GABAa:
IPSPsTA().IPSPmeanUni -- fit double exponential
GABAb?
AMPA:
fit excitatory analogue of IPSPsTA().IPSPmeanUni
NMDA?
Connectivity? (adjacency matrices and kernels)
IPSP/Field relations: PhaseSyncTA()
Spike/Field relations: SpikeFieldCohTA()

------------------------------------------------------------
Distributions (in cells, different LFPs, in layers, ...) -- noted by Tallie
Cell
    Model parameters (parameter distributions)
    Ih size
    Ih time course
    Resistance
    Behavior (distributions to match)
    Spike width
    Spike amp
    ISI @ thresh
    Instant Freq on step
    AHP category (e.g., shape, size, length, ...)
    STO freq at thresh
    STO at high drive

Network
    Model parameters (parameter distributions)
    IPSP/EPSP width & amp (double exponential fit)
    Behavior (distributions to match)
    IPSP/EPSP rhythmicity
    IPSP/EPSP Phase
%}




