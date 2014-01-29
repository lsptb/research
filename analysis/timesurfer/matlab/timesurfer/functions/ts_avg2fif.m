function ts_avg2fif(avg_data,template_fif,outstem,evcode_flag,forceflag)
%function ts_avg2fif(avg_data,template_fif,outstem,evcode_flag,forceflag)
%
% saves each condition in avg_data structure as fif file
%
%  Required parameters:
%   avg_data: timesurfer's avg_data structure
%   template_fif: full pathname of fif file to be used as template
%   outstem: file stem for output files
%
%  Optional parameters:
%   evcode_flag: [0|1] whether to append filenames with event codes (1)
%     or condition numbers (0)
%    {default: 0}
%   forceflag : [0|1] whether to overwrite existing fif files
%    {default: 1}
%
% created:  11/24/05 by Don Hagler
% last mod: 03/26/09 by Jason Sherfey
% last mod: 08/05/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,3)) return; end;
if ~exist('evcode_flag','var') || isempty(evcode_flag), evcode_flag=0; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag=1; end;

[fpath,fstem]=fileparts(outstem);
if ~isempty(fpath) && ~exist(fpath,'dir')
  [s,m] = mkdir(fpath);
  if ~s
    error('failed to create output directory %s:\n%s',fpath,m);
  end;
end;

% get info from online average
if ~exist(template_fif,'file')
  error('template fif file %s not found',template_fif);
end;

megmodel('head',[0,0,0],template_fif);
[na,ki,nu,ct,t]=channames(template_fif);

% save each condition as a fif file
nconds = length(avg_data.averages);
for j=1:nconds
  if evcode_flag
    evcode = avg_data.averages(j).event_code;
    outfile = sprintf('%s_event%d.fif',outstem,evcode);
  else
    outfile = sprintf('%s_cond%d.fif',outstem,j);
  end;
  if ~exist(outfile,'file') || forceflag
    B = double(avg_data.averages(j).data);
    sfreq = avg_data.sfreq;
    t0 = avg_data.averages(j).time(1);
    % added by JSS on 26-Mar-2009
    if size(B,1) < length(na)
      [sel1 sel2] = match_str({avg_data.sensor_info.label},na);
      na = na(sel2);
    end
    savefif(outfile,B,na,sfreq,t0);
  end;
end

