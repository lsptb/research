function parms = getCharParms(source,type,varargin)
% purpose: get parameters for cell characterization functions
% inputs:
% source: 
%   if sim: model spec structure
%   if exp: parmfile
% type: 'sim' or 'exp'
% see also: CharHyperpolStepTA, CharDepolStepTA, CharDepolTonicSpikesTA
% created on 13-Feb-2014 by JSS

cfg = mmil_args2parms( varargin, ...
                         {  'mechanism','iStepProtocol',[],... 
                            'cellfld','entities',[],...
                            'cellind',1,[],...
                            'prestep_time_offset',-.01,[],...
                         }, false);

parms = [];
switch type
  case 'sim' % extract parameters from model specification
    if isstruct(source)
      parms = parsespec(source,cfg);
    else
      error('parm source for simulated data must be the model specification structure.');
    end
  case 'exp' % extract parameters from excel spreadsheet
    if ischar(source) && exist(source,'file')
      parms = parsexls(source,cfg);
    else
      error('parm source for experimental data must be an excel filename');
    end
  otherwise
    error('unknown source type; must be ''sim'' or ''exp''.');
end

function parms = parsespec(spec,cfg)
s = spec.(cfg.cellfld)(cfg.cellind);
p = s.mechs(strcmp(cfg.mechanism,s.mechanisms)).params;

tonic_injected_current = 0;
offset_voltage = 0;
baseline_start_sec = p.bltime/2/1000;%.05;
baseline_length_sec = p.bltime/2/1000;%.05;
prestep_time_offset = cfg.prestep_time_offset; % start section before pulse onset

parms.CharHyperpolStepTA.offset_voltage = offset_voltage;
parms.CharHyperpolStepTA.tonic_injected_current = tonic_injected_current;
parms.CharHyperpolStepTA.sections_label_num = 1:p.nsections;
parms.CharHyperpolStepTA.sections_start_sec = prestep_time_offset + baseline_start_sec + baseline_length_sec;
parms.CharHyperpolStepTA.sections_length_sec = p.isi*p.nsteps/1000;
parms.CharHyperpolStepTA.baseline_start_sec = baseline_start_sec;
parms.CharHyperpolStepTA.baseline_length_sec = baseline_length_sec;

parms.CharDepolStepTA.offset_voltage = offset_voltage;
parms.CharDepolStepTA.tonic_injected_current = tonic_injected_current;
parms.CharDepolStepTA.sections_label_num = 1:p.nsections;
parms.CharDepolStepTA.sections_start_sec = prestep_time_offset + p.nsections*p.isi*p.nsteps/1000 + baseline_start_sec + baseline_length_sec;
parms.CharDepolStepTA.sections_length_sec = p.isi*p.nsteps/1000;
parms.CharDepolStepTA.baseline_start_sec = baseline_start_sec;
parms.CharDepolStepTA.baseline_length_sec = baseline_length_sec;
parms.CharDepolStepTA.step_length = p.steptime/1000;

parms.CharDepolTonicSpikesTA.bpFiltParms = [5 80];
parms.CharDepolTonicSpikesTA.Notch = [];
parms.CharDepolTonicSpikesTA.PlotsYN = 'n';
parms.CharDepolTonicSpikesTA.offset_voltage = offset_voltage;
parms.CharDepolTonicSpikesTA.tonic_injected_current = tonic_injected_current;
parms.CharDepolTonicSpikesTA.sections_label_num = 1;
parms.CharDepolTonicSpikesTA.sections_start_sec = 2*(p.nsections*p.isi*p.nsteps/1000) + baseline_start_sec + baseline_length_sec;
parms.CharDepolTonicSpikesTA.sections_length_sec = p.tonictime/1000;

function parms = parsexls(file,cfg)
% use Tallie's function for loading parameters
parms = GetParmsTA(file);


