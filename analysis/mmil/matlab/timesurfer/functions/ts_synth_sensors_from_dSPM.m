function synth_data = ts_synth_sensors_from_dSPM(stc,prefix,time)
%function synth_data = ts_synth_sensors_from_dSPM(stc,[prefix],[time])
%
% Required Input:
%   stc: source time courses (nsources x time points)
%
% Optional Inupt:
%   prefix: prefix of ts_dSPM output files
%     {default: 'dSPM'}
%   time: time vector (should be in seconds)
%     if empty or ommitted, will use time vector in avg_data struct
%     saved by ts_dSPM
%     {default: []}
%
% Output:
%   synth_data: synthesized data in avg_data structure
%
% Created:  08/01/07 by Don Hagler
% Last Mod: 08/05/09 by Don Hagler
%

synth_data = [];

if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist('prefix','var') | isempty(prefix), prefix = 'dSPM'; end;
if ~exist('time','var'), time = []; end;

if ~isempty(time)
  time = rowvec(time);
  if size(time,2)~=size(stc,2)
    error('length of time vector (%d) does not match time points in stc',...
      size(time,2),size(stc,2));
  end;
end;

matfile=sprintf('matfiles/%s_parms.mat',prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

matfile=sprintf('matfiles/%s_inverse.mat',prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

matfile=sprintf('matfiles/%s_avg_data.mat',prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

G_norm = ts_gain_xyz2norm(G_xyz,...
  parms.lh_dip_info,parms.rh_dip_info,...
  parms.lh_dec_dips,parms.rh_dec_dips,...
  parms.trans);

if size(G_norm,2) ~= size(stc,1)
  error('number of sources in stc (%d) does not match G_norm (%d)\n',...
    size(stc,1),size(G_norm,2));
end;

tmp_data = G_norm*stc;

synth_data = avg_data;
synth_data.averages = synth_data.averages(1);
synth_data.averages.data = zeros(size(synth_data.averages.data,1),...
  size(stc,2));
synth_data.averages.stdev = zeros(size(synth_data.averages.data,1),...
  size(stc,2));
if isempty(time)
  synth_data.averages.time = synth_data.averages.time(1:size(stc,2));
else
  synth_data.averages.time = time;
end;
synth_data.averages.data(parms.goodchans,:) = tmp_data;

