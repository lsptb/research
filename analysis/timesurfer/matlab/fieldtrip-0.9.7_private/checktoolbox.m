function [hastoolbox] = checktoolbox(toolbox, add_to_path);

% CHECKTOOLBOX tests whether an external toolbox is installed
% Optionally it will try to determine the path to the toolbox and
% install it automatically.
% 
% Use as
%   [hastoolbox] = checktoolbox(toolbox, add_to_path);

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: checktoolbox.m,v $
% Revision 1.4  2006/02/07 20:01:39  roboos
% aded biosig and meg-pd (neuromag)
%
% Revision 1.3  2006/01/17 14:05:54  roboos
% added GLNA64 for mentat000
%
% Revision 1.2  2006/01/06 11:39:23  roboos
% added copyrigth and cvs logging, changed some comments
%

% this points the user to the website where he/she can download the toolbox
url = {
  'AFNI'   'http://afni.nimh.nih.gov'
  'DSS'    'http://www.cis.hut.fi/projects/dss'
  'EEGLAB' 'http://www.sccn.ucsd.edu/eeglab'
  'NWAY'   'http://www.models.kvl.dk/source/nwaytoolbox'
  'SPM2'   'http://www.fil.ion.ucl.ac.uk/spm'
  'MEG-PD' 'http://www.kolumbus.fi/kuutela/programs/meg-pd'
  'BIOSIG' 'http://biosig.sourceforge.net'
};

if nargin<2
  % default is not to add the path automatically
  add_to_path = 0;
end

% determine whether the toolbox is installed
toolbox = upper(toolbox);
switch toolbox
case 'AFNI'
  hastoolbox = (exist('BrikLoad') && exist('BrikInfo'));
case 'DSS'
  hastoolbox = exist('dss', 'file') && exist('dss_create_state', 'file');
case 'EEGLAB'
  hastoolbox = exist('runica', 'file');
case 'NWAY'
  hastoolbox = exist('parafac', 'file');
case 'SPM2'
  hastoolbox = (exist('spm_vol') && exist('spm_write_vol') && exist('spm_normalise'));
case 'MEG-PD'
  hastoolbox = (exist('rawdata') && exist('channames'));
case 'BIOSIG'
  hastoolbox = (exist('sopen') && exist('sread'));
otherwise
  warning(sprintf('cannot determine whether the %s toolbox is present', toolbox));
  hastoolbox = 0;
end

% it should be a boolean value
hastoolbox = (hastoolbox~=0);

if ~hastoolbox && add_to_path 
  % try to determine the path of the toolbox
  if strcmp(computer, 'GLNX86') || strcmp(computer, 'GLNXA64')
    % for linux computers in the F.C. Donders Centre
    prefix = '/home/common/matlab/';
  elseif strcmp(computer, 'PCWIN')
    % for windows computers in the F.C. Donders Centre
    prefix = 'h:\common\matlab\';
  else
    prefix = [];
  end
  toolboxpath = [prefix lower(toolbox)];
  if exist(toolboxpath)
    % add the toolbox to the path automatically
    warning(sprintf('adding %s toolbox to your Matlab path', toolbox));
    addpath(toolboxpath);
    hastoolbox = 1;
  else
    % the toolbox is not on the path and cannot be added
    sel = find(strcmp(url(:,1), toolbox));
    if ~isempty(sel)
      msg = sprintf('the %s toolbox is not installed, see %s', toolbox, url{sel, 2});
    else
      msg = sprintf('the %s toolbox is not installed', toolbox);
    end
    error(msg);
  end
end

