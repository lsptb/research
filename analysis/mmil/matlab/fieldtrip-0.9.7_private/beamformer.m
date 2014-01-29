function [dipout] = beamformer(dip, grad, vol, dat, cov, varargin)

% BEAMFORMER scans on pre-defined dipole locations with a single dipole
% and returns the beamformer spatial filter output for a dipole on every
% location.  Dipole locations that are outside the head will return a
% NaN value.
%
% Use as
%   [dipout] = beamformer(dipin, grad, vol, dat, cov, varargin)
% where
%   dipin       is the input dipole model
%   grad        is the gradiometer definition
%   vol         is the volume conductor definition
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   dipout      is the resulting dipole model with all details
%
% The input dipole model consists of
%   dipin.pos   positions for dipole, e.g. regular grid
%   dipin.mom   dipole orientation (optional)
%
% Optional arguments are specified in pairs of a property name and a
% property value
%  'method'          scanning method (see below)
%  'Cr'              cross spectral density between all data channels and the external reference
%  'refdip'          location of dipole with which coherence is computed
%  'refdip'          dipole location(s) that should be suppressed
%  'lambda'          regularisation parameter
%  'powmethod'       can be 'trace' or 'lambda1'
%  'feedback'        give progress indication, can be 'text', 'gui' or 'none' (default)
%  'projectnoise'    project noise estimate through filter,         can be 'yes' or 'no'
%  'keepfilter'      remember the beamformer filter,                can be 'yes' or 'no'
%  'keepleadfield'   remember the forward computation,              can be 'yes' or 'no'
%  'keepmom'         remember the estimated dipole moment,          can be 'yes' or 'no'
%  'keepcsd'         remember the estimated cross-spectral density, can be 'yes' or 'no'
%  'reducerank'      reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%  'normalize'       normalize the leadfield
%
% The scanning method can be
%   'dics'      direct imaging of coherent sources (Gross et al. 2001)
%   'lcmv'      linear constrained minimum variance (van Veen et al. 1997)
%   'pcc'       partial canonical correlation/coherence (Oostenveld & Schoffelen, unpublished)
%
% If the dipole definition only specifies the dipole location, a rotating
% dipole (regional source) is assumed on each location. If a dipole moment
% is specified, its orientation will be used and only the strength will
% be fitted to the data.

% Copyright (C) 2003-2004, Robert Oostenveld
%
% $Log: beamformer.m,v $
% Revision 1.19  2006/03/07 13:11:03  roboos
% added comment, the dics_refdip method will generally not work in combination with a subspace projection
%
% Revision 1.18  2006/01/24 19:58:01  roboos
% renamed submethods for dics from dics_cmc/dics_ccc to dics_refchan/disc_refdip (for consistency with sourceanalysis)
%
% Revision 1.17  2005/11/16 09:56:17  roboos
% added a semicolon to prevent output on screen
%
% Revision 1.16  2005/11/16 09:05:54  roboos
% added optional input argument 'powmethod' (can be 'lambda1' or 'trace'), which determines how the power is computed from the source covariance or csd matrix
%
% Revision 1.15  2005/11/08 11:03:54  roboos
% implemented support for normalize and reducerank using compute_leadfield
% switched from optarg structure to keyval function
% removed obsolete supdip code, since that did not work satisfactorily
%
% Revision 1.14  2005/10/25 08:55:24  roboos
% added some feedback for subspace projection
% changed indentation and some whitespace
%
% Revision 1.13  2005/06/17 08:51:19  roboos
% moved implementation of keepmom and keepcsd from external function into beamformer for more memory efficiency
%
% Revision 1.12  2005/03/31 14:21:35  roboos
% changed noise estimation, take teh maximum of the smallest singular value and lambda (i.e. lambda is the noise floor in case of rank deficiencies)
% shuffled the order in the handling of the common preparation (e.g. noise estimation and rank detection) at the beginning of lcmv and dics
%
% Revision 1.11  2005/02/18 12:50:37  roboos
% fixed bug in computation of leadfield for refdip (due to name change ref->optarg.refdip)
%
% Revision 1.10  2005/02/16 15:30:54  roboos
% fixed bug in the output pos (contained only the inside pos and not the outside)
%
% Revision 1.9  2005/02/16 15:24:21  roboos
% added support for method=pcc using beamformer_pcc function
% replaced the handling of the variable input arguments by an optarg struct
%
% Revision 1.8  2005/02/09 10:54:00  roboos
% previous change in eps was wrong, it needed to be scaled with max(s) as well. Also changed the tolerance to 10 times instead of 2 times the default.
%
% Revision 1.7  2005/02/09 10:45:43  roboos
% anged eps(s) into eps, for compatibility with matlab 61 and 65. This assumes
% that the computations are done in double precision
%
% Revision 1.6  2005/02/08 12:01:16  roboos
% added own implementation of pinv, same as default Matlab but with different tolerance.
% This hopefully solves the numerical inaccuracies that we have encountered on matlab 6.5 and 7.0
%
% Revision 1.5  2004/10/27 16:16:37  roboos
% renamed grid.lbex matrix for subspace projection into grid.subspace (in correspondence with precompute_leadfield)
% transposed the subspace projection matrix
% renamed all occurences of "lbex" (as part of variable names) into "subspace"
%
% Revision 1.4  2004/10/21 17:39:46  roboos
% removed subfunction show_progress, now makes use of external progress
% function which supports more ways of giving feedback on the progress
% and which also handles gui updates in a smarter fashion
%
% Revision 1.3  2004/09/28 09:16:57  roboos
% fixed bug for lcmv+lbex in using the size of the original covariance for regularization
%
% Revision 1.2  2004/09/23 15:04:48  roboos
% iimplemented support for lbex projection in lcmv and dics
% cleaned up spaces and tabs
%
% Revision 1.1  2004/09/14 09:02:48  roboos
% renamed meg_dipole_scan into beamformer, functionally equivalent to revision 1.33
%
% Revision 1.33  2004/08/23 09:31:45  roboos
% changed layout, added output support for source covariance (although not implemented yet), changed conj(x') into ctranspose(x)

if mod(nargin-5,2)
  % the first 5 arguments are fixed, the other arguments should come in pairs
  error('invalid number of optional arguments');
end

% these optional settings do not have defaults
Cr            = keyval('Cr',            varargin);
Pr            = keyval('Pr',            varargin);
refdip        = keyval('refdip',        varargin);
reducerank    = keyval('reducerank',    varargin);
normalize     = keyval('normalize',     varargin);
powmethod     = keyval('powmethod',     varargin); % the default for this is set below
% these optional settings have defaults
method        = keyval('method',        varargin); if isempty(method  ),      error('no method specified'); end
lambda        = keyval('lambda',        varargin); if isempty(lambda  ),      lambda = 0;                   end
feedback      = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end
projectnoise  = keyval('projectnoise',  varargin); if isempty(projectnoise),  projectnoise = 'yes';         end
keepfilter    = keyval('keepfilter',    varargin); if isempty(keepfilter),    keepfilter = 'no';            end
keepleadfield = keyval('keepleadfield', varargin); if isempty(keepleadfield), keepleadfield = 'no';         end
keepmom       = keyval('keepmom',       varargin); if isempty(keepmom),       keepmom = 'yes';              end
keepcsd       = keyval('keepcsd',       varargin); if isempty(keepcsd),       keepcsd = 'no';               end
% convert the yes/no arguments to the corresponding logical values
projectnoise  = strcmp(projectnoise,  'yes');
keepfilter    = strcmp(keepfilter,    'yes');
keepleadfield = strcmp(keepleadfield, 'yes');
keepcsd       = strcmp(keepcsd,       'yes');
keepmom       = strcmp(keepmom,       'yes');

% set the default method by which power is computed from the source covariance/csd
if isempty(powmethod)
  if strcmp(method, 'dics')
    % use the largest singular value of the covariance/csd matrix
    powmethod = 'lambda1';
  elseif strcmp(method, 'lcmv')
    % use the trace of the covariance/csd matrix
    powmethod = 'trace';
  end
end

% use these two logical flags instead of doing the string comparisons each time again
powtrace   = strcmp(powmethod, 'trace');
powlambda1 = strcmp(powmethod, 'lambda1');

% rename this input variable for convenience
if strcmp(method, 'dics')
  Cf = cov;
else
  Cy = cov;
end

% The optional dipole moment can be formatted either as a Nx3 matrix where each row
% specifies a single dipole, or as a vector of length Nx3, where the first 3 values
% correspond with dipole 1, the second three with dipole 2, ...
if isfield(dip, 'mom')
  if ~all(size(dip.mom)==size(dip.pos))
    reshape(dip.mom, 3, size(dip.pos,1))';
  end
end

if ~isempty(Cr)
  % ensure that the cross-spectral density with the reference signal is a column matrix
  Cr = Cr(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside') & ~isfield(dip, 'outside');
  [dip.inside, dip.outside] = find_inside_vol(dip.pos, vol);
elseif isfield(dip, 'inside') & ~isfield(dip, 'outside');
  dip.outside    = setdiff(1:ndip, dip.inside);
elseif ~isfield(dip, 'inside') & isfield(dip, 'outside');
  dip.inside     = setdiff(1:ndip, dip.outside);
end

% select only the dipole positions inside the brain for scanning
dip.posorig = dip.pos;
dip.pos = dip.pos(dip.inside, :);
if isfield(dip, 'mom')
  dip.mom = dip.mom(dip.inside, :);
end
if isfield(dip, 'leadfield')
  dip.leadfield = dip.leadfield(dip.inside);
end
if isfield(dip, 'filter')
  dip.filter = dip.filter(dip.inside);
end
if isfield(dip, 'subspace')
  dip.subspace = dip.subspace(dip.inside);
end

% count the number of dipole positions inside the brain
ndip = size(dip.pos,1);
ntime = size(dat,2);
nchan = size(dat,1);

progress('init', feedback, 'scanning grid');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start the scanning with the proper metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(method)

  case {'partcanoncorr', 'pcc'}
    % this is experimental code that is implemented in a separate function
    dipout = beamformer_pcc(dip, grad, vol, dat, cov, varargin{:});

  case 'dics'
    % dics has the following sub-methods, which depend on the function input arguments
    % power only, cortico-muscular coherence and cortico-cortical coherence
    if ~isempty(Cr) & ~isempty(Pr) & isempty(refdip)
      % compute cortico-muscular coherence, using reference cross spectral density
      submethod = 'dics_refchan';
    elseif isempty(Cr) & isempty(Pr) & ~isempty(refdip)
      % compute cortio-cortical coherence with a dipole at the reference position
      submethod = 'dics_refdip';
    elseif isempty(Cr) & isempty(Pr) & isempty(refdip)
      % only compute power of a dipole at the grid positions
      submethod = 'dics_power';
    else
      error('invalid combination of input arguments for dics');
    end

    % do the preparation for dics that is common to all submethods
    if isfield(dip, 'mom')
      error('dics not yet supported for fixed dipole models');
    end
    isrankdeficient = (rank(Cf)<size(Cf,1));
    if isrankdeficient & ~isfield(dip, 'filter')
      warning('cross-spectral density matrix is rank deficient')
    end
    if projectnoise
      % estimate the noise power, which is further assumed to be equal and uncorrelated over channels
      if isrankdeficient
        % estimated noise floor is equal to or higher than lambda
        noise = lambda;
      else
        % estimate the noise level in the covariance matrix by the smallest singular value
        noise = svd(Cf);
        noise = noise(end);
        % estimated noise floor is equal to or higher than lambda
        noise = max(noise, lambda);
      end
    end
    % the inverse only has to be computed once for all dipoles
    invCf = pinv(Cf + lambda * eye(size(Cf)));
    if isfield(dip, 'subspace')
      fprintf('using subspace projection\n');
      % remember the original data prior to the voxel dependant subspace projection
      Cf_pre_subspace = Cf;
      if strcmp(submethod, 'dics_refchan')
        Cr_pre_subspace = Cr;
        Pr_pre_subspace = Pr;
      end
    end

    switch submethod
      case 'dics_refchan'  % as submethod
        % compute cortico-muscular coherence, using reference cross spectral density
        for i=1:size(dip.pos,1)
          if isfield(dip, 'leadfield')
            % reuse the leadfield that was previously computed
            lf = dip.leadfield{i};
          else
            % compute the leadfield
            lf = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
          end
          if isfield(dip, 'subspace')
            % do subspace projection of the forward model
            lf = dip.subspace{i} * lf;
            % the cross-spectral density becomes voxel dependent due to the projection
            Cr    = dip.subspace{i} * Cr_pre_subspace;
            Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
            invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf))) * dip.subspace{i}');
          end
          if isfield(dip, 'filter')
            % use the provided filter
            filt = dip.filter{i};
          else
            % construct the spatial filter, use PINV/SVD to cover rank deficient leadfield
            filt = pinv(lf' * invCf * lf) * lf' * invCf;
          end
          % compute the source parameters for a dipole on the current position
          if powlambda1
            pow = lambda1(filt * Cf * ctranspose(filt));                   % Gross eqn. 4, 5 and 8
          elseif powtrace
            pow = real(trace(filt * Cf * ctranspose(filt)));
          end
          csd = filt*Cr;                                                   % Gross eqn. 6
          if powlambda1
            coh = lambda1(csd)^2 / (pow * Pr);                             % Gross eqn. 9
          elseif powtrace
            coh = real(trace(csd))^2 / (pow * Pr);
          end
          dipout.pow(i) = pow;
          dipout.coh(i) = coh;
          if keepcsd
            dipout.csd{i} = csd;
          end
          if projectnoise
            if uselambda1
              dipout.noise(i) = noise * lambda1(filt * ctranspose(filt));
            elseif usetrace
              dipout.noise(i) = noise * real(trace(filt * ctranspose(filt)));
            end
          end
          if keepfilter
            dipout.filter{i} = filt;
          end
          if keepleadfield
            dipout.leadfield{i} = lf;
          end
          progress(i/ndip);
        end

      case 'dics_refdip'  % as submethod
        % compute cortio-cortical coherence with a dipole at the reference position
        lf1 = compute_leadfield(refdip, grad, vol, 'reducerank', reducerank, 'normalize', normalize);
        % construct the spatial filter for the first (reference) dipole location, using the non-suppressed coss spectral density
        % use PINV/SVD to cover rank deficient leadfield
        filt1 = pinv(lf1' * invCf * lf1) * lf1' * invCf;
        % compute the power at the first dipole location, using the non-suppressed coss spectral density
        if powlambda1
          Pref = lambda1(filt1 * Cf * ctranspose(filt1));
        elseif powtrace
          Pref = real(trace(filt1 * Cf * ctranspose(filt1)));
        end
        for i=1:size(dip.pos,1)
          if isfield(dip, 'leadfield')
            % reuse the leadfield that was previously computed
            lf2 = dip.leadfield{i};
          else
            % compute the leadfield
            lf2 = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
          end
          if isfield(dip, 'subspace')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % FIXME: the subspace projection is only done on the scandip (lf2), but not on the refdip (lf1),
            % therefore the rotated forward model differs between chandip and refdip and the subspace projection
            % is likely to fail or give nonsense results
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % do subspace projection of the forward model
            lf = dip.subspace{i} * lf;
            % the cross-spectral density becomes voxel dependent due to the projection
            Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
            invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf))) * dip.subspace{i}');
          end
          % construct the spatial filter for the second dipole location
          % use PINV/SVD to cover rank deficient leadfield
          filt2 = pinv(lf2' * invCf * lf2) * lf2' * invCf;
          csd = filt1 * Cf * ctranspose(filt2);                % Gross eqn. 4, compute the cross spectral density between the two dipoles
          if powlambda1
            pow = lambda1(filt2 * Cf * ctranspose(filt2));     % compute the power at the second dipole location
          elseif powtrace
            pow = real(trace(filt2 * Cf * ctranspose(filt2))); % compute the power at the second dipole location)
          end
          if powlambda1
            coh = lambda1(csd)^2 / (pow * Pref);               % compute the coherence between the first and second dipole
          elseif powtrace
            coh = real(trace((csd)))^2 / (pow * Pref);         % compute the coherence between the first and second dipole
          end
          dipout.pow(i) = pow;
          dipout.coh(i) = coh;
          if keepcsd
            dipout.csd{i} = csd;
          end
          if projectnoise
            % project noise on dipole location that is scanned, not the reference location
            if powlambda1
              dipout.noise(i) = noise * lambda1(filt2 * ctranspose(filt2));
            elseif powtrace
              dipout.noise(i) = noise * real(trace(filt2 * ctranspose(filt2)));
            end
          end
          if keepleadfield
            dipout.leadfield{i} = lf2;
          end
          progress(i/ndip);
        end

      case 'dics_power'  % as submethod
        % only compute power of a dipole at the grid positions
        for i=1:size(dip.pos,1)
          if isfield(dip, 'leadfield')
            % reuse the leadfield that was previously computed
            lf = dip.leadfield{i};
          else
            % compute the leadfield
            lf = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
          end
          if isfield(dip, 'subspace')
            % do subspace projection of the forward model
            lf = dip.subspace{i} * lf;
            % the cross-spectral density becomes voxel dependent due to the projection
            Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
            invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf))) * dip.subspace{i}');
          end
          if isfield(dip, 'filter')
            % use the provided filter
            filt = dip.filter{i};
          else
            % construct the spatial filter, use PINV/SVD to cover rank deficient leadfield
            filt = pinv(lf' * invCf * lf) * lf' * invCf;              % Gross eqn. 3
          end
          csd = filt * Cf * ctranspose(filt);                         % Gross eqn. 4 and 5
          if powlambda1
            dipout.pow(i) = lambda1(csd);                             % Gross eqn. 8
          elseif powtrace
            dipout.pow(i) = real(trace(csd));
          end
          if keepcsd
            dipout.csd{i} = csd;
          end
          if projectnoise
            dipout.noise(i) = noise * lambda1(filt * ctranspose(filt));
          end
          if keepfilter
            dipout.filter{i} = filt;
          end
          if keepleadfield
            dipout.leadfield{i} = lf;
          end
          progress(i/ndip);
        end
    end % switch submethod

  case 'lcmv'
    if isfield(dip, 'mom')
      error('lcmv not yet supported for fixed dipole models');
    end
    isrankdeficient = (rank(Cy)<size(Cy,1));
    if isrankdeficient & ~isfield(dip, 'filter')
      warning('covariance matrix is rank deficient')
    end
    if projectnoise
      % estimate the noise power, which is further assumed to be equal and uncorrelated over channels
      if isrankdeficient
        % estimated noise floor is equal to or higher than lambda
        noise = lambda;
      else
        % estimate the noise level in the covariance matrix by the smallest singular value
        noise = svd(Cy);
        noise = noise(end);
        % estimated noise floor is equal to or higher than lambda
        noise = max(noise, lambda);
      end
    end
    % the inverse only has to be computed once for all dipoles
    invCy = pinv(Cy + lambda * eye(size(Cy)));
    if isfield(dip, 'subspace')
      fprintf('using subspace projection\n');
      % remember the original data prior to the voxel dependant subspace projection
      dat_pre_subspace = dat;
      Cy_pre_subspace  = Cy;
    end

    for i=1:size(dip.pos,1)
      if isfield(dip, 'leadfield')
        % reuse the leadfield that was previously computed
        lf = dip.leadfield{i};
      else
        % compute the leadfield
        lf = compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
      end
      if isfield(dip, 'subspace')
        % do subspace projection of the forward model
        lf = dip.subspace{i} * lf;
        % the data and the covariance become voxel dependent due to the projection
        dat   = dip.subspace{i} * dat_pre_subspace;
        Cy    = dip.subspace{i} * (Cy_pre_subspace + lambda * eye(size(Cy_pre_subspace))) * dip.subspace{i}';
        invCy = pinv(dip.subspace{i} * (Cy_pre_subspace + lambda * eye(size(Cy_pre_subspace))) * dip.subspace{i}');
      end
      if powlambda1
        dipout.pow(i) = lambda1(pinv(lf' * invCy * lf));
      elseif powtrace
        dipout.pow(i) = trace(pinv(lf' * invCy * lf));                     % van Veen eqn. 23
      end
      if isfield(dip, 'filter')
        % use the provided filter
        filt = dip.filter{i};
      else
        % construct the spatial filter, use PINV/SVD to cover rank deficient leadfield
        filt = pinv(lf' * invCy * lf) * lf' * invCy;
      end
      % estimate the instantaneous dipole moment at the current position
      if keepmom & ~isempty(dat)
        dipout.mom{i} = filt * dat;
      end
      if projectnoise
        if powlambda1
          dipout.noise(i) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          dipout.noise(i) = noise * trace(filt * ctranspose(filt));
        end
      end
      if keepfilter
        dipout.filter{i} = filt;
      end
      if keepleadfield
        dipout.leadfield{i} = lf;
      end
      progress(i/ndip);
    end

  otherwise
    error(sprintf('unknown method %s', method));
end

progress('close');

dipout.inside  = dip.inside;
dipout.outside = dip.outside;
dipout.pos     = dip.posorig;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'leadfield')
  dipout.leadfield(dipout.inside)  = dipout.leadfield;
  dipout.leadfield(dipout.outside) = {nan};
end
if isfield(dipout, 'filter')
  dipout.filter(dipout.inside)  = dipout.filter;
  dipout.filter(dipout.outside) = {nan};
end
if isfield(dipout, 'mom')
  dipout.mom(dipout.inside)  = dipout.mom;
  dipout.mom(dipout.outside) = {nan};
end
if isfield(dipout, 'pow')
  dipout.pow(dipout.inside)  = dipout.pow;
  dipout.pow(dipout.outside) = nan;
end
if isfield(dipout, 'noise')
  dipout.noise(dipout.inside)  = dipout.noise;
  dipout.noise(dipout.outside) = nan;
end
if isfield(dipout, 'csd')
  dipout.csd(dipout.inside)  = dipout.csd;
  dipout.csd(dipout.outside) = {nan};
end
if isfield(dipout, 'noisecsd')
  dipout.noisecsd(dipout.inside)  = dipout.noisecsd;
  dipout.noisecsd(dipout.outside) = {nan};
end
if isfield(dipout, 'cov')
  dipout.cov(dipout.inside)  = dipout.cov;
  dipout.cov(dipout.outside) = {nan};
end
if isfield(dipout, 'coh')
  dipout.coh(dipout.inside)  = dipout.coh;
  dipout.coh(dipout.outside) = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to obtain the largest singular value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = lambda1(x);
% determine the largest singular value, which corresponds to the power along the dominant direction
s = svd(x);
s = s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the
% standard Matlab function, except that the default tolerance is twice as
% high.
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.19 $  $Date: 2006/03/07 13:11:03 $
%   default tolerance increased by factor 2 (Robert Oostenveld, 7 Feb 2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = pinv(A,varargin)
[m,n] = size(A);
if n > m
  X = pinv(A',varargin{:})';
else
  [U,S,V] = svd(A,0);
  if m > 1, s = diag(S);
  elseif m == 1, s = S(1);
  else s = 0;
  end
  if nargin == 2
    tol = varargin{1};
  else
    tol = 10 * max(m,n) * max(s) * eps;
  end
  r = sum(s > tol);
  if (r == 0)
    X = zeros(size(A'),class(A));
  else
    s = diag(ones(r,1)./s(1:r));
    X = V(:,1:r)*s*U(:,1:r)';
  end
end

